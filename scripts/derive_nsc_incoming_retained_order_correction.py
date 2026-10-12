#!/usr/bin/env python3
"""Record/replay an explicit retained group12 source correction."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_retained_order_correction import prepare_retained_correction,owned_middle_source
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets

OUTPUT=ROOT/'results/development/nsc-incoming-retained-order-correction.json'
SOURCES=('src/recursive_horizons/nsc_incoming_retained_order_correction.py',
         'scripts/derive_nsc_incoming_retained_order_correction.py',
         'tests/test_nsc_incoming_retained_order_correction.py',
         'docs/nsc-incoming-retained-order-correction.md')
INPUTS=('results/development/nsc-incoming-spectral-source.json',
        'results/development/nsc-incoming-middle-bound.json',
        'results/development/nsc-mode-resolved-cauchy-state.json',
        'src/recursive_horizons/nsc_incoming_middle_order_correction.py',
        'src/recursive_horizons/nsc_incoming_source_tail.py',
        'src/recursive_horizons/nsc_incoming_state_moments.py',
        'src/recursive_horizons/nsc_incoming_middle_bound.py',
        'src/recursive_horizons/nsc_incoming_joint_constraints.py',
        'src/recursive_horizons/nsc_incoming_cauchy_jets.py')


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def signature():return{p:sha(p) for p in(*SOURCES,*INPUTS)}


def prepare():
    before=signature();arrays,metadata=prepare_retained_correction(ROOT)
    if before!=signature():raise ValueError('correction dependency changed during preparation')
    metadata['signature']=before
    arrays['metadata_json']=np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-incoming-retained-order-correction.{digest}.npz'
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content-address collision')
    if not path.exists():path.write_bytes(raw)
    return path


def replay(path):
    owned=owned_middle_source(ROOT)
    with np.load(path,allow_pickle=False) as f:arrays={k:f[k].copy() for k in f.files}
    metadata=json.loads(arrays['metadata_json'].tobytes())
    if metadata['signature']!=signature():raise ValueError('correction producer changed')
    for item in metadata['input_payloads']:
        if sha(item['path'])!=item['sha256']:raise ValueError('upstream source payload changed')
    factor=metadata['factor_per_sign'];corrections={};replay_error=0.
    for key,row in metadata['reports'].items():
        if '/' not in key:continue
        value=factor*np.einsum('n,nv->v',arrays[key+'/weights'],arrays[key+'/kernels'])
        replay_error=max(replay_error,float(np.max(abs(value-row['correction']))))
        corrections[key]=value
    change=sum(corrections[f'{s}/refined48/dps70'] for s in(1,-1))
    baseline=np.zeros(4);thermal=np.zeros(4)
    for sign in(1,-1):
        prefix=f'{sign}/original16/'
        baseline+=factor*np.einsum('n,nv->v',arrays[prefix+'weights'],arrays[prefix+'kernels'])
        thermal+=factor*np.einsum('n,nv->v',arrays[prefix+'weights'],arrays[prefix+'thermal_kernels'])
    maxima={'artifact_contraction':replay_error,
            'quadrature':max(metadata['quadrature_indicator_not_bound']),
            'precision':max(metadata['precision_indicator_not_bound']),
            'prefix16':max(row['coefficient_prefix_residual'] for k,row in metadata['reports'].items() if '/' in k),
            'source16_convention':max(metadata['reports'][str(s)]['source16_vacuum_convention'] for s in(1,-1)),
            'source16_decomposition':max(metadata['reports'][str(s)]['source16_decomposition'] for s in(1,-1))}
    tolerances={'artifact_contraction':3e-24,'quadrature':3e-20,'precision':3e-35,
                'prefix16':3e-40,'source16_convention':3e-11,'source16_decomposition':3e-11}
    failures={k:v for k,v in maxima.items() if v>tolerances[k]}
    if failures:raise ArithmeticError(str(failures))
    return {'schema':'NSC-INCOMING-GROUP12-ORDER-CORRECTION-v1','accountable_author':'Douglas Ek',
            'status':'PASS: explicit group12 numerical-order correction; physical remainder separate',
            'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
            'payload':{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()},
            'input_payloads':metadata['input_payloads'],'group':12,'interval':metadata['interval'],
            'old_order':16,'new_order':24,'kernel_order':['rho','p_parallel','T01','p_perp'],
            'old16_middle_approximant':baseline.tolist(),'vacuum_order24_minus16':change.tolist(),
            'corrected_middle_approximant':(baseline+change).tolist(),
            'retained_thermal_insertion':thermal.tolist(),
            'action_gradient_correction':source_action_gradient(incoming_cauchy_jets(),change).tolist(),
            'maxima':maxima,'tolerances':tolerances,'failures':failures,
            'thermal_scattering_bound':metadata['thermal_scattering_bound'],
            'physical_order24_error_bound':None,'physical_error_status':'separate order24 defect certificate required',
            'scope':{'physical_state_action_scales_changed':False,'adiabatic_reference_changed':False,
                     'original_source_replaced':False,'metric_evolution':False,'constraints_solved':False,
                     'old_mode_generators_rerun':False,'numerical_indicators_are_rigorous_bounds':False},
            'reproducer':'python3 scripts/derive_nsc_incoming_retained_order_correction.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.prepare:
        if OUTPUT.exists():raise FileExistsError('existing correction record is not overwritten')
        path=prepare()
    else:
        old=json.loads(OUTPUT.read_text());path=ROOT/old['payload']['path']
        if sha(old['payload']['path'])!=old['payload']['sha256']:raise ValueError('correction artifact changed')
    data=replay(path)
    if args.prepare:OUTPUT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(old,data)
    print(json.dumps({k:data[k] for k in('status','vacuum_order24_minus16','action_gradient_correction','maxima')},indent=2))


if __name__=='__main__':main()
