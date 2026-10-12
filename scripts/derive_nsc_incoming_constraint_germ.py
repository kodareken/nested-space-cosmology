#!/usr/bin/env python3
"""Record structural constraint relation and new independently checked responses."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_constraint_germ import structural_constraint_plane

OUTPUT=ROOT/'results/development/nsc-incoming-constraint-germ.json'
SOURCES=('src/recursive_horizons/nsc_incoming_constraint_germ.py','scripts/derive_nsc_incoming_constraint_germ.py',
         'tests/test_nsc_incoming_constraint_germ.py','docs/nsc-incoming-constraint-germ.md')
INPUTS=('results/development/nsc-incoming-joint-constraints.json',
        'results/development/nsc-incoming-local-constraints.json',
        'results/development/nsc-reference-band-bulk.json',
        'src/recursive_horizons/nsc_incoming_reference_response.py',
        'src/recursive_horizons/nsc_incoming_local_constraints.py',
        'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
        'src/recursive_horizons/nsc_light_restoration_action.py',
        'src/recursive_horizons/nsc_spherical_local_history.py',
        'src/recursive_horizons/nsc_spatial_reference_symbol.py',
        'src/recursive_horizons/nsc_reference_band_action.py')


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def signature():return{p:sha(p) for p in(*SOURCES,*INPUTS)}


def authenticate():
    for p in INPUTS[:3]:
        data=read(p)
        for key in('source_hashes','input_hashes'):
            for path,digest in data.get(key,{}).items():
                if sha(path)!=digest:raise ValueError('response input changed: '+path)
        if 'payload' in data and sha(data['payload']['path'])!=data['payload']['sha256']:
            raise ValueError('response input artifact changed')


def prepare(directory):
    authenticate()
    controls={'u':json.loads((directory/'nsc-incoming-independent-lapse-response.json').read_text()),
              'mixed':json.loads((directory/'nsc-incoming-mixed-germ-response.json').read_text())}
    for control in controls.values():
        for p,digest in control['producer_hashes'].items():
            if sha(p)!=digest:raise ValueError('response producer changed: '+p)
    raw=(json.dumps({'signature':signature(),'controls':controls},indent=2,sort_keys=True)+'\n').encode()
    digest=hashlib.sha256(raw).hexdigest();path=ROOT/f'results/development/artifacts/nsc-incoming-constraint-germ.{digest}.json'
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content-address collision')
    if not path.exists():path.write_bytes(raw)
    return path


def reference_replay(reference):
    total=np.zeros((4,4));first_order=0.
    for sample in reference['quadrature_samples']:
        insertion=np.asarray(sample['insertions']);weights=np.asarray(sample['weights']);k=np.asarray(sample['momenta'])
        if not np.array_equal(k[:len(k)//2],-k[len(k)//2:]):raise ValueError('paired full momentum line required')
        first_order=max(first_order,float(np.max(abs(insertion[0]))))
        total+=sample['multiplicity']/(2*np.pi)*np.einsum('n,onb->ob',weights,insertion)
    result=total.sum(axis=0)
    error=float(np.max(abs(result-reference['reference_action_gradient_change'])))
    return result,error,first_order


def replay(path):
    authenticate();payload=json.loads(path.read_text())
    if payload['signature']!=signature():raise ValueError('constraint-germ dependency changed')
    controls=payload['controls'];u=controls['u'];mixed=controls['mixed']
    old=read(INPUTS[0]);baseline=read(INPUTS[1])['baseline'];base=np.asarray(baseline['local_action_gradient'])
    expectations={'u':{('r',3,0):.01},'mixed':{('r',3,0):.02,('r',1,1):-.005}}
    slot_roundoff=0.;actual_changes={}
    for name,control in controls.items():
        actual={(row['field'],row['normal_order'],row['spatial_order']):row['delta'] for row in control['probe']}
        if actual.keys()!=expectations[name].keys():raise ValueError('owned response plane changed')
        slot_roundoff=max(slot_roundoff,max(abs(actual[k]-expectations[name][k]) for k in actual))
        if slot_roundoff>3e-15:raise ValueError('owned response values changed beyond derivative-slot rounding')
        actual_changes[name]=actual
    u16,err16,first16=reference_replay(u['references']['16'])
    u24,err24,first24=reference_replay(u['references']['24'])
    mixref,errmix,firstmix=reference_replay(mixed['reference'])
    response_u=np.asarray(u['local']['local_action_gradient'])-base+u24[:2]
    response_v=np.asarray(old['local_probe']['local_action_gradient'])-base+np.asarray(old['reference_controls']['n24']['reference_action_gradient_change'])[:2]
    response_mixed=np.asarray(mixed['local']['local_action_gradient'])-base+mixref[:2]
    probe_u=actual_changes['u']['r',3,0]
    mixed_u=actual_changes['mixed']['r',3,0];mixed_v=actual_changes['mixed']['r',1,1]
    cu=response_u[0]/probe_u;cv2=response_v[0]/(.01**2);cv=response_v[1]/.01
    predicted=np.array([cu*mixed_u+cv2*mixed_v**2,cv*mixed_v])
    difference=response_mixed-predicted
    coefficients={'c_u':float(cu),'c_v2':float(cv2),'c_v':float(cv)}
    residuals={'normal_slot_roundoff':slot_roundoff,'reference_node_replay':max(err16,err24,errmix),
               'reference16_24_difference':float(np.max(abs(u16-u24))),
               'homogeneous_shift_zero_identity':float(abs(response_u[1])),
               'independent_mixed_control':float(np.max(abs(difference))),
               'first_reference_order_unchanged':max(first16,first24,firstmix),
               'local_resolution_indicator':max(u['local']['maximum_numerical_indicator'],mixed['local']['maximum_numerical_indicator']),
               'local_Euler_identity':max(u['local']['Euler_bulk_identity_residual'],mixed['local']['Euler_bulk_identity_residual'])}
    if max(residuals.values())>3e-11:raise ArithmeticError('independent response verification failed')
    return {'schema':'NSC-INCOMING-CONSTRAINT-GERM-v1','accountable_author':'Douglas Ek',
            'status':'PASS: structural plane and numerically independent responses; physical root OPEN',
            'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
            'payload':{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()},
            'structure':structural_constraint_plane(),'numerical_coefficients':coefficients,
            'baseline_Jacobian_from_structure':[[float(cu),0.],[0.,float(cv)]],
            'numerical_Jacobian_determinant':float(cu*cv),'rank_evidence':'two numerically resolved independent action responses; not an interval nonsingularity certificate',
            'u_probe_response':response_u.tolist(),'v_probe_response':response_v.tolist(),
            'mixed_control_response':response_mixed.tolist(),'mixed_control_prediction':predicted.tolist(),
            'mixed_control_residual':difference.tolist(),'residuals':residuals,'numerical_tolerance':3e-11,
            'local_probe_records':{'u':u['local'],'mixed':mixed['local']},
            'reference_integrability':{'fixed_grade0_and_grade1_metric_data':True,'Delta_P0_and_P1':'0',
                                       'first_changed_order':2,'reference_difference':'O(abs(k)^-3)',
                                       'raw_vertices':'O(abs(k))','insertion_tail':'O(abs(k)^-2)',
                                       'inventory':'same finite retained33 groups; exact LLL fixed nonzero-k charts',
                                       'owner':'results/development/nsc-reference-band-bulk.json'},
            'scope':{'matter_source_cancels_in_response_difference':True,'physical_state_or_couplings_refit':False,
                     'physical_initial_data_selected':False,'numerical_root_selected':False,
                     'source_accuracy_certified':False,'Cauchy_surface_constraints_solved':False,
                     'parent_preparation_or_endpoint_completed':False,'extended_stationarity':'OPEN',
                     'metric_evolution':False,'Gamma_rest_assigned':False},
            'reproducer':'python3 scripts/derive_nsc_incoming_constraint_germ.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--check',action='store_true');parser.add_argument('--controls',type=Path,default=Path('/tmp'));args=parser.parse_args()
    if args.prepare:
        if OUTPUT.exists():raise FileExistsError('existing constraint-germ record is not overwritten')
        path=prepare(args.controls)
    else:
        old=read(OUTPUT.relative_to(ROOT));path=ROOT/old['payload']['path']
        if sha(old['payload']['path'])!=old['payload']['sha256']:raise ValueError('response artifact changed')
    result=replay(path)
    if args.prepare:OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(old,result)
    print(json.dumps({k:result[k] for k in('status','numerical_coefficients','numerical_Jacobian_determinant','residuals')},indent=2))


if __name__=='__main__':main()
