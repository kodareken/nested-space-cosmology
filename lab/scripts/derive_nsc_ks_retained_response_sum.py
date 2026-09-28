#!/usr/bin/env python3
"""Resumable finite retained-source sum for the already declared nonzero control."""
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
import derive_nsc_ks_retained_response_control as P
from recursive_horizons.nsc_ks_energy_propagator import evolve_energy_propagator,KSEnergyPropagator
from recursive_horizons.nsc_ks_batched_constraints import KSConstraintAccumulator,KSUpstreamBatch,bind_upstream_batch,signed_incoming_pair
from recursive_horizons.nsc_ks_source_inventory import ReferenceSourcePanel
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_evolved_incoming_constraints import compatible_history_slots,surface_geometry_response
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

DIRECTORY=ROOT/'results/development/artifacts/nsc-ks-retained-response-sum'
OUTPUT=ROOT/'results/development/nsc-ks-retained-response-sum.json'
OWNED=('scripts/derive_nsc_ks_retained_response_sum.py','docs/nsc-ks-retained-response-sum.md',
       'src/recursive_horizons/nsc_ks_energy_propagator.py',*P.OWNED)
OP_KEYS=('energies','z','reference','difference','difference_z','tangent','tangent_z')
FAMILY_CAP=90.


def digest(path):return P.digest(path)


def context():
    # The operator/source controls are authenticated and reused; no old PDE run.
    op,meta,_,coeff,baseline,hashes=P.inputs()
    if digest(P.PAYLOAD)!=json.loads(P.OUTPUT.read_text())['payload']['sha256']:
        raise ValueError('pilot source response payload changed')
    inventory=json.loads((ROOT/P.INVENTORY).read_text())
    with np.load(ROOT/inventory['payload']['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
    grouped={}
    for name,p in meta['panels'].items():
        if p['angular_magnitude']==0:continue
        key=(p['group'],p['angular_sign'])
        grouped.setdefault(key,[]).append(ReferenceSourcePanel(name,*key,p['mass'],p['angular_magnitude'],
            *(a[name+'/'+k] for k in ('energies','weights','covariance','amplitudes_at_one')),p['provenance']))
    if len(grouped)!=60:raise ValueError('sixty nonzero-angular positive-source families required')
    signature={**hashes,**{p:digest(p) for p in dict.fromkeys(OWNED)},
        str(P.OUTPUT.relative_to(ROOT)):digest(P.OUTPUT),str(P.PAYLOAD.relative_to(ROOT)):digest(P.PAYLOAD)}
    return op,meta,coeff,baseline,grouped,signature


def paths(key):
    stem=f'family-{key[0]:02d}-{key[1]:+d}'
    return DIRECTORY/(stem+'.json'),DIRECTORY/(stem+'.npz')


def summands(operator,panels,meta,coeff,saved=None):
    family=P.O.C.R.C.P.family(P.O.C.R.C.P.ALPHA)
    acc=KSConstraintAccumulator(family,np.zeros(2),coeff,(operator.z[0],operator.z[-1]),local_error_budget())
    arrays={};receipts=[]
    group,sign=panels[0].group,panels[0].angular_sign
    for panel in panels:
        for part in panel.batches(P.BATCH_ROWS):
            prefix='upstream/'+part.name
            if saved is None:
                batch=bind_upstream_batch(part,rho_up=operator.rho_up,**P.O.C.R.TIGHTER)
                receipt={k:getattr(batch,k) for k in ('panel_name','original_panel','rows','group',
                    'angular_sign','energy_sign','mass','angular','rho_up','preparation_digest')}
                receipt.update(continuation=dict(batch.continuation),provenance=dict(batch.provenance))
            else:
                receipt=saved[1][len(receipts)]
                src=FixedSourcePreparation.from_signed_blocks(part.energies,part.weights,part.covariance)
                batch=KSUpstreamBatch(**receipt,source=src,initial_columns=saved[0][prefix])
            receipts.append(receipt);arrays[prefix]=batch.initial_columns
            state=operator.apply(batch.source,batch.initial_columns)
            negative=part.negative_partner(meta['config'])
            negative_source=FixedSourcePreparation.from_signed_blocks(negative.energies,negative.weights,negative.covariance)
            pos,pbatch,neg,nbatch=signed_incoming_pair(state,batch,negative_source)
            acc.add(pos,pbatch,{**meta['channels'][group],'angular_sign':sign})
            acc.add(neg,nbatch,{**meta['channels'][group],'angular_sign':-sign})
    assembled=acc.finalize()
    arrays['paired_matter_change']=sum(assembled['family_corrections'].values())
    arrays['paired_matter_tangent']=sum(assembled['family_matter_tangents'].values())
    for name,value in assembled['family_corrections'].items():arrays['correction/'+name]=value
    for name,value in assembled['family_matter_tangents'].items():arrays['tangent/'+name]=value
    return arrays,receipts


def new_family(key,op_base,meta,coeff,panels,signature):
    target_json,target_npz=paths(key)
    if target_json.exists() or target_npz.exists():raise FileExistsError('family output already exists')
    cpu=time.process_time()
    def stop(*_):raise TimeoutError('one source family exceeded90 CPU seconds')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,FAMILY_CAP)
    try:
        if key==(14,1):
            # Reuse the complete actual-source pilot; no repeated preparation.
            with np.load(P.PAYLOAD,allow_pickle=False) as f:saved={k:f[k] for k in f.files}
            pmeta=json.loads(saved.pop('metadata_json').tobytes())
            arrays={k:v for k,v in saved.items() if k.startswith(('upstream/','correction/','tangent/'))}
            arrays['paired_matter_change']=sum(v for k,v in arrays.items() if k.startswith('correction/'))
            arrays['paired_matter_tangent']=sum(v for k,v in arrays.items() if k.startswith('tangent/'))
            receipts=pmeta['upstream_batches'];operator=op_base;new_solves=0
        else:
            family=P.O.C.R.C.P.family(P.O.C.R.C.P.ALPHA)
            first=panels[0]
            endpoint=320. if key[0] in (10,11,12,31,32) else 160.
            operator=evolve_energy_propagator((0.,endpoint),48,family,
                op_base.binding.computational_z,op_base.z,first.mass,first.angular,op_base.rho_up,
                axial_support=op_base.binding.axial_support,**P.O.C.R.TIGHTER)
            arrays,receipts=summands(operator,panels,meta,coeff);new_solves=1
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    for k in OP_KEYS:arrays['operator/'+k]=np.asarray(getattr(operator,k))
    arrays['metadata_json']=np.frombuffer(json.dumps({'upstream_batches':receipts},sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays)
    record={'schema':'NSC-KS-RETAINED-RESPONSE-FAMILY-v1','positive_family':list(key),
        'operator_energy_interval':operator.interval,'operator_digest':operator.digest,
        'source_rows':sum(len(p.energies) for p in panels),'source_panels':[p.name for p in panels],
        'quadrature_measure':sum(float(p.weights.sum()) for p in panels),
        'signed_families':sorted(k.split('/',1)[1] for k in arrays if k.startswith('correction/')),
        'matter_change_maxima':np.max(abs(arrays['paired_matter_change']),axis=0).tolist(),
        'CPU_seconds':time.process_time()-cpu,'family_CPU_cap':FAMILY_CAP,'new_operator_solves':new_solves,
        'upstream_batches':len(receipts),'source_hashes':signature,
        'physical_local_gate':'OPEN','source_and_numerical_error_bound':None,
        'payload':{'path':str(target_npz.relative_to(ROOT)),'sha256':sha256(raw).hexdigest(),'bytes':len(raw)}}
    target_npz.write_bytes(raw);target_json.write_text(json.dumps(record,sort_keys=True,indent=2,allow_nan=False)+'\n')
    return record


def read_family(key,op_base,meta,coeff,panels,signature,replay=False):
    path,_=paths(key);r=json.loads(path.read_text())
    if r['source_hashes']!=signature or r['positive_family']!=list(key):raise ValueError('family signature/identity changed')
    if r['source_rows']!=sum(len(p.energies) for p in panels) or r['source_panels']!=[p.name for p in panels]:
        raise ValueError('family source coverage changed')
    expected_measure=320. if key[0] in (10,11,12,31,32) else 160.
    if abs(r['quadrature_measure']-expected_measure)>1e-11:raise ValueError('family quadrature inventory incomplete')
    if digest(r['payload']['path'])!=r['payload']['sha256']:raise ValueError('family payload changed')
    with np.load(ROOT/r['payload']['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
    first=panels[0]
    operator=KSEnergyPropagator(tuple(r['operator_energy_interval']),*(a['operator/'+k] for k in OP_KEYS),
        first.mass,first.angular,op_base.rho_up,op_base.binding,{})
    if operator.digest!=r['operator_digest']:raise ValueError('family operator binding changed')
    if sorted(k.split('/',1)[1] for k in a if k.startswith('correction/'))!=r['signed_families']:
        raise ValueError('signed family coverage changed')
    if not np.array_equal(sum(v for k,v in a.items() if k.startswith('correction/')),a['paired_matter_change']):
        raise ValueError('signed response sum changed')
    if not np.array_equal(sum(v for k,v in a.items() if k.startswith('tangent/')),a['paired_matter_tangent']):
        raise ValueError('signed tangent sum changed')
    if replay:
        pm=json.loads(a['metadata_json'].tobytes())
        current,_=summands(operator,panels,meta,coeff,(a,pm['upstream_batches']))
        for k,v in current.items():
            if not np.array_equal(v,a[k]):raise ValueError('family replay changed: '+k)
    return r,a['paired_matter_change'],a['paired_matter_tangent']


def assemble(op,meta,coeff,baseline,grouped,signature,replay=False):
    records=[];dG=np.zeros((len(op.z),2));dJ=np.zeros((1,len(op.z),2));signed=[]
    for key in sorted(grouped):
        if not paths(key)[0].exists():continue
        record,delta,tangent=read_family(key,op,meta,coeff,grouped[key],signature,replay)
        records.append(record);dG+=delta;dJ+=tangent;signed.extend(record['signed_families'])
    if len(signed)!=len(set(signed)):raise ValueError('signed family counted twice')
    family=P.O.C.R.C.P.family(P.O.C.R.C.P.ALPHA)
    slots,dslots=compatible_history_slots(family,op.z,1);geometry=surface_geometry_response(slots,dslots,coeff)
    residual=np.asarray(baseline['baseline']['action_gradient_approximant'])+geometry['action_gradient_change']+dG
    full=len(records)==len(grouped)
    record={'schema':'NSC-KS-RETAINED-RESPONSE-SUM-v1','accountable_author':'Douglas Ek',
        'status':'OPEN: '+('complete retained finite-source control; uncertified physical gate' if full else 'partial retained finite-source control'),
        'completed_positive_families':len(records),'required_positive_families':len(grouped),
        'all_retained_nonzero_angular_families_included':full,
        'positive_energy_rows':sum(r['source_rows'] for r in records),
        'signed_family_ids':sorted(signed),'z':op.z.tolist(),
        'summed_matter_change':dG.tolist(),'summed_matter_tangent':dJ.tolist(),
        'raw_constraint_approximant':residual.tolist(),
        'raw_history_jacobian':(geometry['action_gradient_tangent']+dJ).tolist(),
        'constraint_maxima':np.max(abs(residual),axis=0).tolist(),
        'matter_change_maxima':np.max(abs(dG),axis=0).tolist(),
        'family_records':[{'path':str(paths(tuple(r['positive_family']))[0].relative_to(ROOT)),
            'sha256':digest(paths(tuple(r['positive_family']))[0])} for r in records],
        'CPU_seconds':sum(r['CPU_seconds'] for r in records),'new_operator_solves':sum(r['new_operator_solves'] for r in records),
        'source_hashes':signature,'state_law':'C_Sigma[g]=U_g C_up U_gdagger; original fixed source',
        'scope':{'local_interval':'S(1)+[.12,.18]','physical_local_gate':'OPEN',
            'full_source_accuracy':None,'node_and_interpolation_error':None,'changed_history_UV_bound':None,
            'between_node_error':None,'constraint_root_claimed':False,'metric_timestep':False,
            'baseline_and_geometry_counted_once':True,'energy_folding_factor':1,
            'analytic_zero_radius_response_groups':[0,13,23],'ell0_baseline_retained':True,
            'history_amplitude':P.O.C.R.C.P.ALPHA,'stress_drift_subtracted':False},
        'reproducer':'python3 scripts/derive_nsc_ks_retained_response_sum.py --check'}
    return record


def run(cpu_budget,max_new):
    if cpu_budget<=0 or max_new<=0:raise ValueError('positive work bounds required')
    op,meta,coeff,baseline,grouped,signature=context();DIRECTORY.mkdir(parents=True,exist_ok=True)
    start=time.process_time();new=0
    for key in sorted(grouped):
        if paths(key)[0].exists():
            read_family(key,op,meta,coeff,grouped[key],signature);continue
        if time.process_time()-start>=cpu_budget or new>=max_new:break
        r=new_family(key,op,meta,coeff,grouped[key],signature);new+=1
        print(json.dumps({k:r[k] for k in ('positive_family','source_rows','matter_change_maxima','CPU_seconds')}),flush=True)
        progress=assemble(op,meta,coeff,baseline,grouped,signature)
        (DIRECTORY/'progress.json').write_text(json.dumps(progress,sort_keys=True,indent=2,allow_nan=False)+'\n')
    result=assemble(op,meta,coeff,baseline,grouped,signature)
    if result['all_retained_nonzero_angular_families_included']:
        OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--cpu-budget',type=float,default=1200.);p.add_argument('--max-new',type=int,default=60)
    args=p.parse_args()
    if args.run:result=run(args.cpu_budget,args.max_new)
    else:
        op,meta,coeff,baseline,grouped,signature=context()
        result=assemble(op,meta,coeff,baseline,grouped,signature,replay=True)
        path=OUTPUT if OUTPUT.exists() else DIRECTORY/'progress.json'
        if result!=json.loads(path.read_text()):raise ValueError('source sum replay changed')
    print(json.dumps({k:result[k] for k in ('status','completed_positive_families','required_positive_families',
        'positive_energy_rows','constraint_maxima','matter_change_maxima','CPU_seconds')},indent=2))
