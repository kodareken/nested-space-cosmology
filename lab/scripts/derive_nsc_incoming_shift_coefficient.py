#!/usr/bin/env python3
"""Prepare/replay the closed c_v coefficient; no old probe/reference runs."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_shift_coefficient import (
    weyl_trace_identity,reference_shift_expressions,integrated_reference_expressions,
    local_compact_identities,local_density_controls,closed_shift_coefficient,
)
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets

OUTPUT='results/development/nsc-incoming-shift-coefficient.json'
SOURCES=('src/recursive_horizons/nsc_incoming_shift_coefficient.py',
         'scripts/derive_nsc_incoming_shift_coefficient.py',
         'tests/test_nsc_incoming_shift_coefficient.py','docs/nsc-incoming-shift-coefficient.md')
RECORDS=('results/development/nsc-incoming-joint-constraints.json',
         'results/development/nsc-incoming-local-constraints.json',
         'results/development/nsc-incoming-constraint-germ.json',
         'results/development/nsc-mode-resolved-cauchy-state.json',
         'results/development/nsc-reference-band-bulk.json')
OWNERS=('src/recursive_horizons/nsc_spatial_reference_symbol.py',
        'src/recursive_horizons/nsc_general_ks_reference.py',
        'src/recursive_horizons/nsc_incoming_reference_response.py',
        'src/recursive_horizons/nsc_incoming_cauchy_jets.py',
        'src/recursive_horizons/nsc_incoming_local_constraints.py',
        'src/recursive_horizons/nsc_light_restoration_action.py',
        'src/recursive_horizons/nsc_spherical_local_history.py',
        'src/recursive_horizons/nsc_magnetic_light_reference.py',
        'src/recursive_horizons/nsc_incoming_lapse_coefficient.py')


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def signature():return {p:digest(p) for p in (*SOURCES,*RECORDS,*OWNERS)}
def read(p):return json.loads((ROOT/p).read_text())


def inputs():
    records=[read(p) for p in RECORDS]
    for record in records:
        for key in ('source_hashes','input_hashes'):
            for p,h in record.get(key,{}).items():
                if digest(p)!=h:raise ValueError('shift coefficient input changed: '+p)
        if 'payload' in record and digest(record['payload']['path'])!=record['payload']['sha256']:
            raise ValueError('stored response artifact changed')
    return records


def distance(value,bounds,scale=1.):
    lo,hi=sorted([scale*bounds[0],scale*bounds[1]])
    return max(0.,lo-value,value-hi)


def stored_response_controls(records,closed):
    joint,local,germ,inventory,_=records
    with np.load(ROOT/joint['payload']['path'],allow_pickle=False) as a:
        meta=json.loads(a['metadata_json'].tobytes()); control=meta['controls']['n24']
        changes=meta['probe_changes']
        if len(changes)!=1 or (changes[0]['field'],changes[0]['normal_order'],changes[0]['spatial_order'])!=('r',1,1):
            raise ValueError('stored same-plane v response required')
        v=changes[0]['delta']
        domain=incoming_cauchy_jets();n=domain.normal_geometry()
        scalar=lambda x:float(x.value[0,0,0].real)
        A2=scalar(domain.fields[2].derivative(t=2))/n['intrinsic_a']
        R2=scalar(domain.fields[3].derivative(t=2))/n['intrinsic_r']
        formulas=reference_shift_expressions()
        functions={order:sp.lambdify(formulas['symbols'],expr,'numpy',cse=True) for order,expr in formulas['kernels'].items()}
        node_errors=dict.fromkeys(('1','2','3','4'),0.)
        for index,row in enumerate(control['signed_inventory']):
            k=a[f'n24/{index}/momenta']; actual=a[f'n24/{index}/insertions']
            m=row['mass'];L=row['angular']/n['intrinsic_r'];p=k/n['intrinsic_a'];w=np.sqrt(m*m+L*L+p*p)
            args=(m,L,p,w,n['intrinsic_r'],n['k_parallel'],n['k_perp'],A2,R2)
            for order,f in functions.items():
                predicted=v*np.asarray(f(*args))
                node_errors[str(order)]=max(node_errors[str(order)],float(np.max(abs(predicted-actual[order-1,:,1]))))
    reference={}
    for order in (2,4):
        old=joint['reference_controls']['n24']['formal_order_gradient_changes'][order-1][1]
        reference[str(order)]=distance(old,closed['reference_order_intervals'][str(order)],v)
    channel_errors={}
    base=local['baseline']['channel_action_gradients'];probe=joint['local_probe']['channel_action_gradients']
    for channel in base:
        bounds=closed['local_channel_intervals'].get(channel,[0.,0.])
        channel_errors[channel]=distance(probe[channel][1]-base[channel][1],bounds,v)
    payload=read(germ['payload']['path'])
    mixed_v=next(row['delta'] for row in payload['controls']['mixed']['probe'] if (row['field'],row['normal_order'],row['spatial_order'])==('r',1,1))
    return {'stored_node_residuals':node_errors,'reference_integrated_response_residuals':reference,
        'local_channel_response_residuals':channel_errors,
        'total_v_response_residual':distance(germ['v_probe_response'][1],closed['total_coefficient_interval'],v),
        'independent_mixed_response_residual':distance(germ['mixed_control_response'][1],closed['total_coefficient_interval'],mixed_v),
        'previous_small_probe_coefficient':germ['numerical_coefficients']['c_v'],
        'previous_coefficient_distance':distance(germ['numerical_coefficients']['c_v'],closed['total_coefficient_interval']),
        'stored_probe_v':v,'stored_mixed_v':mixed_v,
        'stored_reference_payload':joint['payload'],'stored_mixed_payload':germ['payload'],
        'new_reference_nodes_prepared':False}


def prepare():
    records=inputs();before=signature();closed=closed_shift_coefficient(records[3]['channels'],records[1]['locked_inputs'])
    ref=integrated_reference_expressions(); raw=reference_shift_expressions()
    proof={'weyl_trace':weyl_trace_identity(),'bloch_normalization_residuals':raw['normalization_residuals'],
        'reference_coefficients':{str(n):str(v) for n,v in ref['compact_coefficients'].items()},
        'reference_moments':ref['moments'],'reference_compact_reduction_residuals':ref['compact_reduction_residuals'],
        'odd_reference_full_line_residual':ref['odd_order_full_line_identity'],
        'local':local_compact_identities()}
    controls=stored_response_controls(records,closed); density=local_density_controls(records[1]['locked_inputs'])
    if signature()!=before:raise ValueError('shift producer/input changed during exact coefficient preparation')
    payload={'schema':'NSC-INCOMING-SHIFT-COEFFICIENT-ARTIFACT-v1','signature':before,
             'proof':proof,'controls':controls,'local_density_controls':density,'directed_coefficient':closed}
    raw=(json.dumps(payload,indent=2,sort_keys=True,allow_nan=False)+'\n').encode();h=sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-incoming-shift-coefficient.{h}.json'
    if path.exists() and path.read_bytes()!=raw:raise ValueError('content-addressed collision')
    if not path.exists():path.write_bytes(raw)
    return path


def make_record(path):
    records=inputs();payload=json.loads(path.read_text())
    if payload['signature']!=signature():raise ValueError('shift coefficient producer or input changed')
    closed=closed_shift_coefficient(records[3]['channels'],records[1]['locked_inputs'])
    if closed!=payload['directed_coefficient']:raise ValueError('directed closed coefficient replay differs')
    controls=payload['controls'];proof=payload['proof']
    zeros=[*proof['bloch_normalization_residuals'],*proof['reference_compact_reduction_residuals'],
           proof['odd_reference_full_line_residual'],*proof['local']['closed_form_residuals'].values(),*proof['local']['plane_residuals'].values()]
    if set(zeros)!={'0'} or not closed['strictly_positive']:raise ArithmeticError('exact coefficient identities or nonzero interval failed')
    residuals={'stored_reference_nodes':max(controls['stored_node_residuals'].values()),
        'stored_reference_order_responses':max(controls['reference_integrated_response_residuals'].values()),
        'stored_local_channel_responses':max(controls['local_channel_response_residuals'].values()),
        'stored_total_v_response':controls['total_v_response_residual'],
        'stored_independent_mixed_response':controls['independent_mixed_response_residual'],
        'original_density_derivatives':payload['local_density_controls']['maximum'],
        'directed_replay':0.}
    failures={k:v for k,v in residuals.items() if not np.isfinite(v) or v>3e-11}
    return {'schema':'NSC-INCOMING-SHIFT-COEFFICIENT-v1','accountable_author':'Douglas Ek',
        'status':'PASS: closed strictly positive c_v; physical initial data and stationarity OPEN' if not failures else 'OPEN: coefficient control mismatch',
        'coefficient':closed,'proof':proof,'controls':controls,'residuals':residuals,'control_tolerance':3e-11,'failures':failures,
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in (*RECORDS,*OWNERS)},
        'payload':{'path':str(path.relative_to(ROOT)),'sha256':sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size},
        'scope':{'coefficient_fixed_by_existing_action':True,'same_full_Weyl_and_2D_local_Euler':True,
            'physical_source_cancels_in_response':True,'source_state_or_scale_refit':False,
            'c_u_proof_rerun':False,'probe_reference_or_mode_generators_run':False,
            'c_v2_newly_certified':False,'physical_initial_data_selected':False,'Cauchy_surface_solution':False,
            'physical_source_accuracy_certified':False,'parent_preparation_or_endpoints_completed':False,
            'extended_stationarity':'OPEN','metric_evolution':False,'Gamma_rest_assigned':False},
        'reproducer':'python3 scripts/derive_nsc_incoming_shift_coefficient.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.prepare:
        result=make_record(prepare());(ROOT/OUTPUT).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:
        old=read(OUTPUT);spec=old['payload']
        if digest(spec['path'])!=spec['sha256']:raise ValueError('coefficient artifact changed')
        result=make_record(ROOT/spec['path'])
        if result!=old:raise ValueError('coefficient record differs from authenticated replay')
    print(json.dumps({'status':result['status'],'interval':result['coefficient']['total_coefficient_interval'],
                      'residuals':result['residuals'],'failures':result['failures']},indent=2))
    return 1 if result['failures'] else 0


if __name__=='__main__':raise SystemExit(main())
