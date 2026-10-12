#!/usr/bin/env python3
"""Convert a completed fine-history field proof to finite-source matter error."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb,acb,acb_mat,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import validate_nsc_ks_fine_history as V
import derive_nsc_ks_fine_trajectory as F
from recursive_horizons import nsc_ks_source_envelope as E
from recursive_horizons.nsc_ks_difference_envelope import KSDifferenceIncoming
from recursive_horizons.nsc_ks_ball_trajectory import BallFourierSegment,complex_ball,exact_upper,restored_upper
from recursive_horizons.nsc_ks_reference_error import spectator_reference_error,homogeneous_reference_norms
from recursive_horizons.nsc_ks_finite_matter_error import source_norm_upper,endpoint_norms,finite_matter_error
from recursive_horizons.nsc_ks_local_constraints import homogeneous_reference_amplitudes,stationary_computational_reference
from recursive_horizons.nsc_ks_matter_difference import source_fixed_matter_difference
from recursive_horizons.nsc_ks_streaming_trajectory import TrajectoryChunkReader

OUTPUT=ROOT/'results/development/nsc-ks-fine-matter-accuracy.json'
OWNED=('scripts/derive_nsc_ks_fine_matter_accuracy.py','docs/nsc-ks-fine-matter-accuracy.md',
       'src/recursive_horizons/nsc_ks_reference_error.py','src/recursive_horizons/nsc_ks_finite_matter_error.py',
       'src/recursive_horizons/nsc_ks_local_constraints.py','src/recursive_horizons/nsc_ks_matter_difference.py',
       'src/recursive_horizons/nsc_evolved_incoming_constraints.py','src/recursive_horizons/nsc_common_subtracted_ks_source.py')
BITS=120
CONTROL_TOLERANCE=1e-11


def digest(path):
    path=Path(path);return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def restore_state(a,meta,family):
    grid=a['computational_z'];w,U=E._sample_axial_profiles(family.directions,grid)
    binding=E.KSEnvelopeBinding(grid,w,U,family.amplitudes,
        tuple(d.inner_radius for d in family.directions),tuple(d.outer_radius for d in family.directions),
        E.usual_axial_support(),meta['rho_up'],1.,**meta['options'])
    if binding.fingerprint!=meta['history_fingerprint']:raise ValueError('fine history binding changed')
    get=lambda key:a['prepared/'+key]
    return KSDifferenceIncoming(a['target_z'],get('columns'),get('column_tangents'),get('axial_columns'),
        get('axial_tangents'),a['source/covariance'],a['source/column_weights'],a['source/energies'],
        meta['parameters']['mass'],meta['parameters']['angular'],a['initial_canonical_columns'],meta['rho_up'],
        binding,meta['preparation'],None,meta['diagnostics'],get('reference_amplitudes'),
        get('envelope_difference'),get('envelope_difference_z'))


def matrix(array):return acb_mat([[complex_ball(v) for v in row] for row in array])


def support_hull(supports):
    if not supports:raise ValueError('at least one axial support required')
    low,high=supports[0]
    for left,right in supports[1:]:
        low=arb.min(low,left);high=arb.max(high,right)
    return low,high


def insertion(Fcol,Fzcol,C,mass,angular,multiplicity):
    D=Fcol*C*Fcol.conjugate().transpose()
    J=(Fzcol*C*Fcol.conjugate().transpose()-Fcol*C*Fzcol.conjugate().transpose())/acb(0,2)
    s1=acb_mat([[0,1],[1,0]]);s2=acb_mat([[0,acb(0,-1)],[acb(0,1),0]])
    s3=acb_mat([[1,0],[0,-1]])
    radius=arb(2).sqrt();axial=(3*arb.pi()/2-4).sqrt()
    return (multiplicity*(mass*(s1*D).trace()-angular/radius*(s2*D).trace()-(s3*J).trace()/axial),
            multiplicity*J.trace())


def node_errors(field,state,reference,reference_z,cache,origin):
    maximum=[arb(0) for _ in range(4)]
    for node,z in enumerate(state.z):
        X=field.value(1,arb(float(z))-origin)
        Xz=field.value(1,arb(float(z))-origin,axial_derivative=1)
        errors=[arb(0) for _ in range(4)]
        for spin in range(2):
            for source,(energy,weight) in enumerate(zip(state.source_energies,state.column_weights)):
                phase=acb(0,-arb(float(energy))*arb(float(z))).exp()
                current=phase*X[spin][source]
                current_z=phase*(Xz[spin][source]-acb(0,arb(float(energy)))*X[spin][source])
                ref=phase*complex_ball(cache[spin,source])*arb(float(weight))
                refz=-acb(0,arb(float(energy)))*ref
                exact=(current,current_z,ref,refz)
                observed=(state.weighted_columns[node,spin,source],state.weighted_axial_columns[node,spin,source],
                          reference[node,spin,source]*weight,reference_z[node,spin,source]*weight)
                for i,(x,y) in enumerate(zip(exact,observed)):
                    errors[i]+=(x-complex_ball(y)).abs_upper()**2
        for i,error in enumerate(errors):maximum[i]=arb.max(maximum[i],error.upper().sqrt().upper())
    return tuple(v.upper() for v in maximum)


def compute():
    if not V.OUTPUT.exists():raise ValueError('all history cells must be validated first')
    validation=json.loads(V.OUTPUT.read_text())
    if not validation['all_segments_covered'] or validation['numerical_field_error_relative_to_saved_upstream'] is None:
        raise ValueError('partial history bounds cannot certify matter accuracy')
    for path,expected in validation['signature'].items():
        if digest(path)!=expected:raise ValueError('validated history code/input changed')
    for item in validation['segment_proofs']:
        if digest(item['path'])!=item['sha256']:raise ValueError('segment proof changed')
    fine=json.loads(F.OUTPUT.read_text())
    if digest(fine['payload']['path'])!=fine['payload']['sha256']:raise ValueError('fine prepared fields changed')
    with np.load(F.FINAL,allow_pickle=False) as f:a={k:f[k] for k in f.files}
    meta=json.loads(a['metadata_json'].tobytes());family=F.T.C.R.C.P.family(F.T.C.R.C.P.ALPHA)
    state=restore_state(a,meta,family)
    reader=TrajectoryChunkReader(F.DIRECTORY);last=reader.load_chunk(len(reader)-1)
    segment=last.segment(len(last)-1)
    cached=homogeneous_reference_amplitudes(state)
    reference,reference_z=stationary_computational_reference(state)
    params=meta['parameters']
    numerical=source_fixed_matter_difference(state,axial_scale=params['axial_scale'],radius=params['radius'],
                                            multiplicity=params['multiplicity'],reference='cached')
    with ctx.workprec(BITS):
        field=BallFourierSegment(segment,state.column_weights,len(a['computational_z']),meta['period_length'],bits=BITS)
        norm,z_norm=endpoint_norms(field,state.source_energies)
        field_error=validation['numerical_field_error_relative_to_saved_upstream']
        epsF=restored_upper(field_error['weighted_F_error_upper']);epsZ=restored_upper(field_error['weighted_F_z_error_upper'])
        size=state.initial_columns.size;Dend=segment.end[size:].reshape(2,len(state.source_energies),len(a['computational_z']))
        supports=[(arb(d.w.center)-arb(d.w.outer),arb(d.w.center)+arb(d.w.outer)) for d in family.directions]
        low,high=support_hull(supports)
        origin=arb(float(a['computational_z'][0]))
        ref_error=spectator_reference_error(state.rho_up,origin,field.length,low,high,Dend[:,:,0],
            state.column_weights,state.source_energies,restored_upper(field_error['row_field_error_upper']),
            reference_pair=(state.reference_amplitudes,cached),bits=BITS)
        eval_errors=node_errors(field,state,reference,reference_z,cached,origin)
        refF,refZ=homogeneous_reference_norms(cached,state.column_weights,state.source_energies)
        gamma=source_norm_upper(state.source_covariance)
        parameters=dict(mass=state.mass,absolute_angular=abs(state.angular),axial_lower=arb(4)/5,
                        radius_lower=arb(7)/5,multiplicity=params['multiplicity'],bits=BITS)
        current=finite_matter_error(norm+eval_errors[0],z_norm+eval_errors[1],epsF+eval_errors[0],epsZ+eval_errors[1],gamma,**parameters)
        ref=finite_matter_error(refF+eval_errors[2],refZ+eval_errors[3],
            restored_upper(ref_error['weighted_reference_F_error_upper'])+eval_errors[2],
            restored_upper(ref_error['weighted_reference_F_z_error_upper'])+eval_errors[3],gamma,**parameters)
        C=matrix(state.source_covariance);arithmetic=[arb(0),arb(0)]
        for node in range(len(state.z)):
            cur=insertion(matrix(state.weighted_columns[node]),matrix(state.weighted_axial_columns[node]),C,
                          arb(state.mass),arb(state.angular),arb(params['multiplicity']))
            old=insertion(matrix(reference[node]*state.column_weights),matrix(reference_z[node]*state.column_weights),C,
                          arb(state.mass),arb(state.angular),arb(params['multiplicity']))
            for j in range(2):
                difference=cur[j]-old[j]-arb(float(numerical['action_gradient_change'][node,j]))
                arithmetic[j]=arb.max(arithmetic[j],difference.abs_upper())
        total=[(restored_upper(current[key])+restored_upper(ref[key])+arithmetic[j]).upper() for j,key in enumerate(('N','beta'))]
        passed=all(v<=arb(CONTROL_TOLERANCE) for v in total)
        result={'schema':'NSC-KS-FINE-MATTER-ACCURACY-v1','accountable_author':'Douglas Ek',
            'status':('PASS' if passed else 'OPEN')+': finite-source numerical matter error; physical gate OPEN',
            'numerical_control_pass':passed,'numerical_control_tolerance':CONTROL_TOLERANCE,
            'current_field_error':current,'reference_field_error':ref,'reference_error_proof':ref_error,
            'node_evaluation_error_upper':[exact_upper(v) for v in eval_errors],
            'arithmetic_and_intrinsic_geometry_error_upper':[exact_upper(v) for v in arithmetic],
            'total_N_beta_error_upper':[exact_upper(v) for v in total],
            'source_matrix_operator_norm_upper':exact_upper(gamma),
            'scope':{'source_family':'14_1','sampled_incoming_nodes':len(state.z),
                'full_source_accuracy':None,'source_energy_UV_bound':None,'constraint_root_claimed':False,
                'physical_local_gate':'OPEN','metric_timestep':False,'stress_drift_subtracted':False},
            'source_hashes':{p:digest(p) for p in OWNED},
            'input_hashes':{str(V.OUTPUT.relative_to(ROOT)):digest(V.OUTPUT),str(F.OUTPUT.relative_to(ROOT)):digest(F.OUTPUT)},
            'reproducer':'python scripts/derive_nsc_ks_fine_matter_accuracy.py --check'}
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args()
    if args.run:
        if OUTPUT.exists():raise FileExistsError('matter accuracy record exists; use --check')
        result=compute();OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:
        saved=json.loads(OUTPUT.read_text());result=compute()
        if result!=saved:raise ValueError('matter accuracy replay differs')
    with ctx.workprec(BITS):display=[float(restored_upper(v)) for v in result['total_N_beta_error_upper']]
    print(json.dumps({'status':result['status'],'N_beta_error_upper_display':display,'tolerance':CONTROL_TOLERANCE},indent=2))
