#!/usr/bin/env python3
"""Reproduce the same-prescription curved Weyl cocycle and its source."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import eigh,expm

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_covariant_operator import (
    smooth_metric,cutoff_response,ultraviolet_subtraction,heat_coefficients)
from recursive_horizons.nsc_influence import canonical_hamiltonian,ground_covariance,influence
from recursive_horizons.nsc_finite_terms import finite_response
from recursive_horizons.nsc_shape_response import StaticAxialMetric
from recursive_horizons.nsc_regulated import RegulatedOperator,OperatorConventions
from recursive_horizons.nsc_weyl_cocycle import (
    conformal_metric,ratio_hamiltonian,ratio_operator,cutoff_functional,
    heat_path_cocycle,operator_variations,directional_derivatives,compensated_lift,
    local_anomaly_cocycle,continuum_identities)
OUTPUT=ROOT/'results/nsc-12-weyl-cocycle.json'


def native(value):
    if isinstance(value,complex):return {'real':float(value.real),'imag':float(value.imag)}
    if isinstance(value,dict):return {k:native(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [native(v) for v in value]
    if isinstance(value,np.ndarray):return native(value.tolist())
    if isinstance(value,np.generic):return value.item()
    return value


def compare_record(expected,actual,path='$'):
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or expected.keys()!=actual.keys():raise RuntimeError('keys differ at '+path)
        for k in expected:compare_record(expected[k],actual[k],path+'/'+k)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(expected)!=len(actual):raise RuntimeError('list differs at '+path)
        for i,(a,b) in enumerate(zip(expected,actual)):compare_record(a,b,f'{path}/{i}')
    elif isinstance(expected,float):
        if isinstance(actual,bool) or not isinstance(actual,(float,int)) or not np.isclose(expected,actual,atol=2e-7,rtol=2e-8):
            raise RuntimeError(f'numeric field differs at {path}: {expected} != {actual}')
    elif type(expected) is not type(actual) or expected!=actual:raise RuntimeError('exact field differs at '+path)


def axial(m):return StaticAxialMetric(m.length,m.lapse,m.radial_scale,m.sphere_radius)


def shift_metric(m,direction,t):
    return replace(m,**{name:getattr(m,name)*np.exp(t*v)
        for name,v in zip(('lapse','radial_scale','sphere_radius'),direction)})


def calculate():
    metric=smooth_metric(16,general=True)
    angle=2*np.pi*metric.x/metric.length
    phi=.04+.07*np.cos(angle);sigma=.03*np.sin(angle)-.02*np.cos(2*angle)
    changed=conformal_metric(metric,phi);twice=conformal_metric(changed,sigma)
    cutoff,kmax,nf=2.,6,48
    raw=[cutoff_functional(m,cutoff,kmax,nf) for m in (metric,changed,twice)]
    direct=raw[1]['energy']-raw[0]['energy']
    path=[{'points':n,'heat_integral':heat_path_cocycle(metric,phi,cutoff,kmax,nf,n)} for n in (4,8)]
    for row in path:row['minus_direct']=row['heat_integral']-direct
    second_path=heat_path_cocycle(changed,sigma,cutoff,kmax,nf,8)
    total_path=heat_path_cocycle(metric,phi+sigma,cutoff,kmax,nf,8)
    composition={'first_path':path[-1]['heat_integral'],'second_path':second_path,
                 'combined_path':total_path,'residual':path[-1]['heat_integral']+second_path-total_path}
    lift=compensated_lift(metric,phi,cutoff,kmax,nf)
    gauge_shift=compensated_lift(conformal_metric(metric,sigma),phi-sigma,cutoff,kmax,nf)
    gauge={'base_energy':lift['base_energy'],'lifted_energy':lift['lifted_energy'],
           'compensating_action':lift['compensating_action'],
           'gauge_shifted_energy':gauge_shift['lifted_energy'],
           'energy_difference':gauge_shift['lifted_energy']-lift['lifted_energy'],
           'Ward_residual':lift['Ward_residual'],
           'phi_gradient_norm':float(np.linalg.norm(lift['phi_gradient'])),
           'metric_gradient_norm':float(np.linalg.norm(lift['metric_gradients'])),
           'lift_is_a_stationary_solution':False,'gauge_interpretation_is_conditional':True}
    h=ratio_hamiltonian(metric);hc=ratio_hamiltonian(changed)
    _,_,cov=ground_covariance(h);_,_,covc=ground_covariance(hc)
    time=.7;u=expm(-1j*time*h);uc=expm(-1j*time*hc)
    history=influence(cov,uc,u)
    state={'Hamiltonian_difference_norm':float(np.linalg.norm(h-hc)),
           'covariance_difference_norm':float(np.linalg.norm(cov-covc)),
           'unitary_difference_norm':float(np.linalg.norm(u-uc)),
           'classical_history_amplitude':history['amplitude'],
           'time':time,'raw_action_change':direct,
           'same_classical_propagation_does_not_imply_measure_invariance':True}
    phase=[]
    for kappa in (1,2):
        d=ratio_operator(metric,.7,kappa);dc=ratio_operator(changed,.7,kappa)
        a=RegulatedOperator(d);b=RegulatedOperator(dc)
        phase.append({'kappa':kappa,'negative_count_before':int(sum(a.values<0)),
                      'negative_count_after':int(sum(b.values<0)),
                      'relative_unwrapped_phase':b.negative_logdet_phase()-a.negative_logdet_phase(),
                      'smallest_absolute_eigenvalue':float(min(min(abs(a.values)),min(abs(b.values))))})
    # Finite full-matrix norms need not converge; compare resolved states and traces.
    grids=[]
    for n in (12,16,24,32):
        m=smooth_metric(n,general=True);s=.04+.07*np.cos(2*np.pi*m.x/m.length)
        m=conformal_metric(m,s);new=ratio_hamiltonian(m);old=canonical_hamiltonian(m)
        e,vec=eigh(new,driver='evd');low=vec[:,abs(e)<3.]
        raw_new=cutoff_functional(m,2.,6,48,False)['energy']
        raw_old=cutoff_response(m,2.,6,48,gradients=False)['energy']
        grids.append({'points':n,'full_H_difference_norm':float(np.linalg.norm(new-old)),
                      'resolved_H_difference_norm':float(np.linalg.norm(low.conj().T@(new-old)@low)),
                      'resolved_energy_window':3.,'resolved_rank':low.shape[1],
                      'new_energy':raw_new,'prior_energy':raw_old,'energy_difference':raw_new-raw_old})
    small=smooth_metric(12,general=True);a=2*np.pi*small.x/small.length
    direction=np.array([.3*np.cos(a),.2*np.sin(a),.4+.1*np.cos(2*a)])
    small_raw=cutoff_functional(small,2.,4,48)
    derivatives=directional_derivatives(small,direction,2.,4,48)
    derivatives['independent_spectral_gradient']=float(np.sum(small_raw['log_gradients']*direction))
    fd=[]
    for step in (.004,.002):
        ep=cutoff_functional(shift_metric(small,direction,step),2.,4,48,False)['energy']
        em=cutoff_functional(shift_metric(small,direction,-step),2.,4,48,False)['energy']
        first=(ep-em)/(2*step);second=(ep-2*small_raw['energy']+em)/step**2
        fd.append({'step':step,'first':first,'second':second,
                   'first_error':first-derivatives['first'],'second_error':second-derivatives['second']})
    local_rows=[]
    for n,scale in ((16,6.),(32,6.),(48,2.),(48,4.),(48,6.),(64,6.),(64,8.)):
        m=smooth_metric(n,general=True);p=.05+.04*np.cos(2*np.pi*m.x/m.length)
        mp=conformal_metric(m,p);angular=int(5*scale*max(mp.sphere_radius));frequencies=96
        e0=cutoff_functional(m,scale,angular,frequencies,False)['energy']
        e1=cutoff_functional(mp,scale,angular,frequencies,False)['energy']
        sub0=ultraviolet_subtraction(m,scale,1.)['energy'];sub1=ultraviolet_subtraction(mp,scale,1.)['energy']
        local=local_anomaly_cocycle(m,p)
        observed=(e1-sub1)-(e0-sub0)
        local_rows.append({'points':n,'cutoff':scale,'angular_max':angular,'frequency_points':frequencies,
                           'phi_formula':'0.05+0.04*cos(2*pi*x/L)',
                           'raw_cocycle':e1-e0,'subtraction_difference':sub1-sub0,
                           'subtracted_cocycle':observed,'local_anomaly_integral':local,
                           'difference_from_anomaly_limit':observed-local,
                           'resolved_spatial_control':n>=48})
    constant=.08
    high_metric=smooth_metric(64,general=True)
    high_phi=.05+.04*np.cos(2*np.pi*high_metric.x/high_metric.length)
    high_changed=conformal_metric(high_metric,high_phi)
    high_angular=int(40*max(high_changed.sphere_radius))
    high_controls=[{'angular_max':high_angular,'frequency_points':96,
                    'subtracted_cocycle':local_rows[-1]['subtracted_cocycle']}]
    for angular,frequencies in ((high_angular,128),(high_angular+8,128)):
        energies=[]
        for m in (high_metric,high_changed):
            energies.append(cutoff_functional(m,8.,angular,frequencies,False)['energy']
                            -ultraviolet_subtraction(m,8.,1.)['energy'])
        high_controls.append({'angular_max':angular,'frequency_points':frequencies,
                              'subtracted_cocycle':energies[1]-energies[0]})
    constant_metric=smooth_metric(32,general=True)
    scaled=cutoff_functional(conformal_metric(constant_metric,constant),2.,20,64,False)['energy']
    raised=cutoff_functional(constant_metric,2.*np.exp(constant),20,64,False)['energy']
    reference=cutoff_functional(constant_metric,2.,20,64,False)['energy']
    heat=heat_coefficients(constant_metric)
    approximation=(heat['a0']*2.**4*np.expm1(4*constant)/4
                   +heat['a2']*2.**2*np.expm1(2*constant)/2+heat['a4']*constant)/(16*np.pi**2)
    rigid={'phi':constant,'points':32,'angular_max':20,'frequency_points':64,'cutoff':2.,
           'scaled_metric_energy':scaled,'raised_cutoff_energy':raised,
           'operator_scaling_identity_residual':scaled-raised,'direct_cocycle':scaled-reference,
           'local_heat_approximation':approximation,'finite_cutoff_remainder':scaled-reference-approximation,
           'proper_time_a0_exponential_denominator':4,'hard_projector_a0_exponential_denominator':8}
    inv_metric=smooth_metric(128,general=True);p=.04+.07*np.cos(2*np.pi*inv_metric.x/inv_metric.length)
    finite0=finite_response(axial(inv_metric));finite1=finite_response(axial(conformal_metric(inv_metric,p)))
    shape=.3*np.cos(2*np.pi*inv_metric.x/inv_metric.length)
    invariant={'C2_energy':float(finite0['energy'][2]),
               'C2_Weyl_difference':float(finite1['energy'][2]-finite0['energy'][2]),
               'C2_radius_shape_response':float(inv_metric.spacing*np.dot(
                   finite0['gradients'][2,2],inv_metric.sphere_radius*shape)),
               'coefficient_selected':False}
    # The independently derived normalization flow must commute with this cocycle.
    from recursive_horizons.nsc_normalization_flow import finite_coefficient_shift
    ratio=2.3;delta=finite_coefficient_shift(ratio)
    f0=finite_response(axial(metric));f1=finite_response(axial(changed))
    uv0=ultraviolet_subtraction(metric,2.,1.);uv1=ultraviolet_subtraction(changed,2.,1.)
    uv0p=ultraviolet_subtraction(metric,2.,ratio);uv1p=ultraviolet_subtraction(changed,2.,ratio)
    c0=direct-(uv1['energy']-uv0['energy']);c1=direct-(uv1p['energy']-uv0p['energy'])
    increment=float(np.dot(delta,f1['energy']-f0['energy']))
    g0=uv0['metric_gradients']-uv0p['metric_gradients']+metric.spacing*np.einsum('k,kij->ij',delta,f0['gradients'])
    g1=uv1['metric_gradients']-uv1p['metric_gradients']+changed.spacing*np.einsum('k,kij->ij',delta,f1['gradients'])
    cocycle_gradient=np.exp(phi)[None,:]*g1-g0
    phi_residual=np.sum(g1*np.array([changed.lapse,changed.radial_scale,changed.sphere_radius]),axis=0)
    combined={'normalization_ratio':ratio,'subtracted_cocycle_before':c0,'subtracted_cocycle_after':c1,
              'derived_finite_increment_difference':increment,'matched_square_residual':c1+increment-c0,
              'metric_gradient_square_max_by_field':np.max(abs(cocycle_gradient),axis=1),
              'phi_gradient_square_max':float(max(abs(phi_residual))),
              'initial_coefficients':None,'ratio_identified_with_Omega':False}
    gates={'ratio_discretization_preserves_Weyl':state['Hamiltonian_difference_norm']<1e-12,
           'same_classical_history':abs(history['amplitude']-1)<2e-12 and state['covariance_difference_norm']<1e-12,
           'covariant_source_local_Ward':max(r['local_Weyl_Ward_residual'] for r in raw)<1e-11,
           'heat_path_matches_action':abs(path[-1]['minus_direct'])<3e-10,
           'heat_path_cocycle_composes':abs(composition['residual'])<3e-10,
           'metric_gradient_and_hessian_agree':abs(fd[-1]['first_error'])<2e-6 and abs(fd[-1]['second_error'])<2e-6,
           'same_continuum_as_prior_scheme':abs(grids[-1]['energy_difference'])<2e-10 and grids[-1]['resolved_H_difference_norm']<1e-10,
           'relative_phase_preserved':all(r['relative_unwrapped_phase']==0 for r in phase),
           'conditional_gauge_lift_verified':abs(gauge['energy_difference'])<1e-11 and gauge['Ward_residual']<1e-11,
           'normalization_and_compensator_square_closes':abs(combined['matched_square_residual'])<1e-12
               and max(combined['metric_gradient_square_max_by_field'])<1e-12 and combined['phi_gradient_square_max']<1e-12,
           'finite_invariant_shape_equation_remains':abs(invariant['C2_Weyl_difference'])<1e-9 and abs(invariant['C2_radius_shape_response'])>1.,
           'local_anomaly_approached':abs(local_rows[-1]['difference_from_anomaly_limit'])<1.3e-5,
           'high_cutoff_quadratures_resolved':max(abs(high_controls[i]['subtracted_cocycle']-high_controls[i-1]['subtracted_cocycle']) for i in (1,2))<1e-8,
           'cutoff6_spatial_refinement':abs(local_rows[-2]['subtracted_cocycle']-local_rows[-3]['subtracted_cocycle'])<1e-8}
    if not all(gates.values()):raise RuntimeError(f'Weyl cocycle gate failed: {gates}')
    sources=('scripts/check_nsc_weyl_cocycle.py','src/recursive_horizons/nsc_weyl_cocycle.py',
             'src/recursive_horizons/nsc_covariant_operator.py','src/recursive_horizons/nsc_finite_terms.py',
             'src/recursive_horizons/nsc_influence.py','src/recursive_horizons/nsc_regulated.py','src/recursive_horizons/nsc_shape_response.py',
             'src/recursive_horizons/nsc_normalization_flow.py','tests/test_nsc_weyl_cocycle.py','docs/nsc-weyl-cocycle.md')
    inputs=('results/nsc-9-covariant-source.json','results/nsc-10-measure-normalization.json',
            'results/nsc-10-influence.json','results/nsc-11-response-matching.json','results/nsc-12-normalization-flow.json')
    return native({'schema':'nsc-weyl-cocycle-v1','artifact_id':'NSC-12-WEYL-COCYCLE',
        'classification':'computed_same_prescription_Weyl_anomaly_integrand_and_conditional_compensator_lift',
        'conventions':{'geometry':'smooth R2 axial circle, a=1, static positive N,q,r, AP symmetric Fourier grids',
          'operator':'H=rho2{N/q,P}/2+rho1*N*kappa/r; D=N^-1/2 beta(omega+iH)N^-1/2',
          'discretization':'ratio form is an explicit finite representation preserving common Weyl covariance; prior records unchanged',
          'regulator':'four-dimensional proper-time modulus only; no unregulated rank integral or separate heat action',
          'units':'hbar=c=a=1; energy and cocycle per coordinate time in inverse a',
          'angular_counting':'representative2N spinor block times4kappa',
          'Weyl':'N,q,r multiply by exp(phi); M,Lambda held fixed, no Omega or zeta identification',
          'conditional_lift':'E_hat[g,phi]=E[g_phi]; phi gauge interpretation conditional, no measure selected'},
        'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'input_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        'development_settings':{'points':16,'cutoff':2.,'angular_max':6,'frequency_points':48,
            'metric_owner':'nsc_covariant_operator.smooth_metric(points=16,general=True)',
            'lapse':metric.lapse,'radial_scale':metric.radial_scale,'sphere_radius':metric.sphere_radius},
        'metric_derivative_settings':{'points':12,'cutoff':2.,'angular_max':4,'frequency_points':48,
            'log_metric_directions_N_q_r':direction},
        'exact_continuum_identities':continuum_identities(),'phi_profile':phi,'sigma_profile':sigma,
        'raw_actions':[r['energy'] for r in raw],'local_Weyl_residuals':[r['local_Weyl_Ward_residual'] for r in raw],
        'direct_cocycle':direct,'independent_heat_path':path,'composition':composition,
        'conditional_lift':gauge,'classical_state_control':state,'relative_phase':phase,
        'discretization_refinement':grids,'metric_derivatives':derivatives,'independent_differences':fd,
        'subtracted_anomaly_limit':local_rows,'high_cutoff_angular_frequency_controls':{
            'points':64,'cutoff':8.,'phi_formula':'0.05+0.04*cos(2*pi*x/L)','rows':high_controls},
        'constant_Weyl_and_profile_weights':rigid,
        'remaining_Weyl_invariant_term':invariant,'normalization_Weyl_square':combined,
        'gate':gates,'nonclaims':{'full_scale_field_measure_derived':False,'fixed_background_average_regularized':False,
            'physical_Weyl_gauge_status_proved':False,'finite_invariant_coefficients_selected':False,
            'full_in_in_compensator_or_relative_scale_action':False,'self_sourced_geometry':False,
            'physical_mass_or_Omega_derived':False,'rigorous_continuum_error_bound':False,
            'new_to_world_Weyl_anomaly_mechanism':False},'terminal':True})


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite record')
    record=calculate()
    if args.check:compare_record(json.loads(OUTPUT.read_text()),record);print('Curved Weyl cocycle, source and normalization square reproduced; all fields checked.')
    elif args.output:
        with args.output.open('x') as f:json.dump(record,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')
        print(args.output)
    else:print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
