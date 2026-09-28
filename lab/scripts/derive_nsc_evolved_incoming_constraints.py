#!/usr/bin/env python3
"""Join saved evolved-state controls to constraints; never propagate fields."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import make_interp_spline

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_evolved_incoming_constraints import (
    evolved_constraint_diagnostic,source_column_matter,column_derivatives)
from recursive_horizons.nsc_evolved_incoming_state import (
    EvolvedIncomingState,FixedSourcePreparation,HistoryBinding,_coverage)
import derive_nsc_evolved_incoming_state as STATE
from derive_nsc_incoming_source_update_v3 import authenticate

OUTPUT='results/development/nsc-evolved-incoming-constraints.json'
SOURCES=('src/recursive_horizons/nsc_evolved_incoming_constraints.py',
         'scripts/derive_nsc_evolved_incoming_constraints.py',
         'tests/test_nsc_evolved_incoming_constraints.py','docs/nsc-evolved-incoming-constraints.md')
RECORDS=('results/development/nsc-evolved-incoming-state.json',
         'results/development/nsc-incoming-source-update-v5.json',
         'results/development/nsc-incoming-surface-coefficients.json',
         'results/development/nsc-incoming-surface-regular-branch.json',
         'results/development/nsc-incoming-fixed-transfer.json')
OWNERS=('src/recursive_horizons/nsc_common_subtracted_ks_source.py',
        'src/recursive_horizons/nsc_prepared_history_jets.py',
        'src/recursive_horizons/nsc_evolved_incoming_state.py',
        'scripts/derive_nsc_evolved_incoming_state.py',
        'scripts/derive_nsc_incoming_source_update_v3.py',
        'results/development/nsc-retarded-compatible-response.json')
MAX=lambda x:float(np.max(abs(x),initial=0.))


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()


def coefficients_from_records(coeff,branch):
    c,b=coeff['coefficients'],branch['branch'];mid=lambda x:float(np.mean(x))
    return {'a':mid(b['intrinsic_intervals']['a']),'r':mid(b['intrinsic_intervals']['r']),
        'Hr':mid(b['intrinsic_intervals']['Hr']),'A0':mid(b['imported_unchanged']['A0']),
        'A1':mid(b['intervals']['A1_total']),'d':mid(b['imported_unchanged']['d']),
        'e':mid(b['imported_unchanged']['e']),'C':mid(c['C_interval']),
        'D':np.r_[0.,[mid(c['polynomial_intervals']['D'][str(j)]) for j in range(1,5)]],
        'F':np.array([mid(c['polynomial_intervals']['F'][str(j)]) for j in range(3)])}


def load_inputs():
    records=[json.loads((ROOT/path).read_text()) for path in RECORDS]
    for record in records:authenticate(record)
    state_record,baseline,coefficient,branch,transfer=records
    if not state_record['status'].startswith('PASS:') or state_record['control']['field_runs']!=3:
        raise ValueError('completed nonzero-history state control required')
    if baseline['partial_error_budget']['full_source_error_bound'] is not None:
        raise ValueError('this control records the current explicit baseline accuracy hole')
    with np.load(ROOT/state_record['payload']['path'],allow_pickle=False) as archive:
        arrays={key:archive[key] for key in archive.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    channel=transfer['actual_source_bindings']['channel']
    if (channel['index']!=14 or channel['compact_mass']!=arrays['source/mass'].item()
            or channel['angular_eigenvalue']!=arrays['source/angular'].item()):
        raise ValueError('same authenticated signed group14_1 operator required')
    source=FixedSourcePreparation(arrays['source/covariance'],arrays['source/column_weights'],arrays['source/energies'])
    return records,arrays,meta,channel,source,coefficients_from_records(coefficient,branch)


def restore_saved_control(arrays,meta,source,name,amplitude):
    """Only called after authentication of the actual evolved-control payload."""
    times,x=arrays['times'],arrays['x'];rho1=int(np.flatnonzero(x==1.)[0])
    binding=HistoryBinding(amplitudes=(amplitude,),inner_radii=(.007,),outer_radii=(.03,),
        times=times,spatial_grid=x,rho1_index=rho1,rho1_metric=arrays['center/rho1_metric'],
        restriction=arrays['center/restriction'],restriction_parameter_derivative=0.,
        S1=meta['support']['S1'],source_digest=source.digest)
    fields=arrays[name+'/columns']
    tangents=arrays['center/column_tangents'][None] if name=='center' else np.zeros((0,*fields.shape),complex)
    return EvolvedIncomingState(arrays['center/z'],fields,tangents,source.covariance,
        source.column_weights,source.energies,binding,_coverage(fields.shape[-1],4),
        {'restored_from_authenticated_evolution':True,'case':name})


def calculate():
    records,a,meta,channel,source,coeff=load_inputs()
    state_record,baseline=records[:2]
    z=a['center/z'];times=a['times']
    # Derivative of a cubic interpolant of the saved trace; raw F values stay
    # unchanged. Continuum derivative error is explicitly not certified.
    Dz=make_interp_spline(z,np.eye(len(z)),k=3).derivative()(z)
    x,phi,data,_=STATE.RETARDED.inputs()
    if not np.array_equal(x,a['x']):raise ValueError('same original upstream grid required')
    ix=int(np.flatnonzero(x==1.)[0]);rows=np.array([ix,len(x)+ix])
    ref0=a['center/restriction']@phi[rows]
    reference=ref0[None]*np.exp(-1j*times[:,None,None]*source.energies[None,None,:])
    reference_z=-1j*source.energies*reference  # stationary computational reference ONLY
    alpha=state_record['control']['amplitude'];h=state_record['control']['finite_difference_step']
    outputs={};states={}
    for name,amplitude in (('center',alpha),('plus',alpha+h),('minus',alpha-h)):
        state=restore_saved_control(a,meta,source,name,amplitude);states[name]=state
        outputs[name]=evolved_constraint_diagnostic(state,STATE.family(amplitude,meta['support']),Dz,
            reference,reference_z,baseline['baseline']['action_gradient_approximant'],coeff,channel,1)
    center=outputs['center'];fd=(outputs['plus']['partial_constraint_diagnostic']-outputs['minus']['partial_constraint_diagnostic'])/(2*h)
    residual=MAX(fd-center['history_derivative_diagnostic'][0])
    state=states['center'];F=state.weighted_columns;dF=state.weighted_column_tangents
    Fz,dFz=column_derivatives(F,dF,Dz)
    params=dict(mass=channel['compact_mass'],angular=channel['angular_eigenvalue'],
                axial_scale=coeff['a'],radius=coeff['r'],multiplicity=channel['copy_count']*channel['degeneracy']/2)
    actual=source_column_matter(F,Fz,source.covariance,dF,dFz,**params)
    old=json.loads((ROOT/'results/development/nsc-retarded-compatible-response.json').read_text())
    authenticate(old)
    with np.load(ROOT/old['payload']['path'],allow_pickle=False) as archive:
        zero_columns=archive['801_64/trace']
        if (not np.array_equal(archive['801_64/times'],times)
                or not np.array_equal(archive['source/source'],a['source/blocks'])
                or not np.array_equal(archive['source/weights'],a['source/weights'])):
            raise ValueError('saved zero-history control must share the same time/source inventory')
    zeros=np.zeros((0,*zero_columns.shape),complex)
    zero_weighted=zero_columns*source.column_weights
    zero_z,_=column_derivatives(zero_weighted,zeros,Dz)
    zero_matter=source_column_matter(zero_weighted,zero_z,source.covariance,zeros,zeros,**params)
    reference_matter=source_column_matter(reference*source.column_weights,reference_z*source.column_weights,
                                         source.covariance,zeros,zeros,**params)
    drift=zero_matter['action_gradient']-reference_matter['action_gradient']
    history_change=actual['action_gradient']-zero_matter['action_gradient']
    drift_split={'reference_field_and_derivative_drift_indicator':np.max(abs(drift),axis=0).tolist(),
        'history_change_at_same_discretization':np.max(abs(history_change),axis=0).tolist(),
        'sum_residual':MAX(drift+history_change-center['evolved_matter_correction']),
        'drift_subtracted_from_main_diagnostic':False,'is_continuum_error_bound':False,
        'saved_zero_history_payload':old['payload']}
    wrong=source_column_matter(F,-1j*source.energies*F,source.covariance,
                               dF,-1j*source.energies*dF,**params)
    diagonal=source_column_matter(F,Fz,np.diag(np.diag(source.covariance)),dF,dFz,**params)
    effects={'dropped_delta_C':MAX(center['evolved_matter_tangent']),
             'outgoing_k_equals_minus_E':MAX(actual['action_gradient_tangent']-wrong['action_gradient_tangent']),
             'dropped_source_coherence':MAX(actual['action_gradient_tangent']-diagonal['action_gradient_tangent'])}
    rejected=[]
    for label,value in (('bare_C0',np.eye(2)/2),('legacy_matter',{'stress_approximant':np.zeros(4)})):
        try:
            evolved_constraint_diagnostic(value,STATE.family(alpha,meta['support']),Dz,
                reference,reference_z,baseline['baseline']['action_gradient_approximant'],coeff,channel,1)
        except TypeError:rejected.append(label)
        else:raise AssertionError('frozen input accepted as evolved state')
    return records,a,meta,center,residual,effects,rejected,MAX(Dz@np.ones(len(z))),MAX(
        np.einsum('ij,jas->ias',Dz,reference)-reference_z),drift_split


def make_record():
    records,arrays,meta,center,residual,effects,rejected,constant_error,ref_derivative_error,drift_split=calculate()
    # Fixed before evaluating a root: a discrete wiring test, NOT 3e-11
    # physical stationarity. Omitted-state controls must be separately visible.
    derivative_tolerance=3e-8
    passed=residual<derivative_tolerance and min(effects.values())>10*residual
    return {'schema':'NSC-EVOLVED-INCOMING-CONSTRAINTS-v1','accountable_author':'Douglas Ek',
        'status':'PASS: evolved-state constraint wiring; full physical constraints OPEN' if passed else
                 'OPEN: integrated derivative or omitted-term discrimination unresolved',
        'state_law':center['state_law'],'control':{'family':'14_1','source_energies':4,
            'amplitude':records[0]['control']['amplitude'],'finite_difference_step':records[0]['control']['finite_difference_step'],
            'field_runs':0,'array_only':True,'axial_derivative':'derivative of cubic interpolant on saved z nodes',
            'source_energy_substitution_used_on_evolved_columns':False,
            'reference':'original stationary canonical columns for fixed computational subtraction only',
            'physical_constraint_tolerance':3e-11,'derivative_control_tolerance':derivative_tolerance},
        'residuals':{'complete_discrete_constraint_derivative':residual,'constant_derivative':constant_error,
                     'reference_axial_derivative_indicator':ref_derivative_error},
        'omitted_term_effects':effects,'frozen_inputs_rejected':rejected,
        'reference_drift_separation':drift_split,
        'diagnostics':{'z':center['z'].tolist(),'partial_constraint_diagnostic':center['partial_constraint_diagnostic'].tolist(),
            'history_derivative_diagnostic':center['history_derivative_diagnostic'].tolist(),
            'evolved_matter_correction':center['evolved_matter_correction'].tolist(),
            'evolved_matter_tangent':center['evolved_matter_tangent'].tolist(),
            'maximum_absolute_partial_gradient':np.max(abs(center['partial_constraint_diagnostic']),axis=0).tolist(),
            'sampled_channel':center['sampled_channel']},
        'error_budget':{'baseline_covered_region_error':records[1]['partial_error_budget']['covered_spectral_regions_action_error_upper'],
            'baseline_full_source_error':None,'changed_history_complete_source_error':None,
            'field_and_axial_derivative_error':None,'coefficient_interval_propagation':None},
        'scope':{'fixed_C0_regression_unchanged':True,'fixed_source_and_evolved_incoming':True,
            'retarded_delta_C_included':True,'old_fixed_C0_roots_relabelled':False,
            'reference_counted_once':True,'reference_state_tangent_subtracted':False,
            'partial_channel_control_only':True,'full_constraint_residual_certified':False,
            'physical_constraint_status':'OPEN','physical_history_selected':False,'metric_timestep':False,
            'new_stress_or_action_term':False,'locked_parameter_refit':False,'publication':False},
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in (*RECORDS,*OWNERS)},
        'reused_artifacts':[records[0]['payload'],drift_split['saved_zero_history_payload']],
        'reproducer':'python3 scripts/derive_nsc_evolved_incoming_constraints.py --check'}


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();record=make_record();path=ROOT/OUTPUT
    if args.write:path.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=record:raise ValueError('evolved constraint receipt differs')
    print(json.dumps({key:record[key] for key in ('status','residuals','omitted_term_effects','error_budget')},indent=2))


if __name__=='__main__':main()
