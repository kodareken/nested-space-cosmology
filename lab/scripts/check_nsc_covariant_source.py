#!/usr/bin/env python3
"""Reproduce the covariant frequency-integrated Dirac metric source."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_covariant_operator import (
    SIGMA1,SIGMA2,SIGMA3,cutoff_response,cutoff_weyl_trace,cylinder_metric,
    cylinder_zero_winding_energy,directional_variation,euclidean_operator,
    frequency_density,heat_coefficients,heat_trace,local_trace_coefficient,
    radial_momentum,smooth_metric,spatial_cutoff_proxy,subtracted_response,
    ultraviolet_subtraction,ultrastatic_reference,ward_summary,
)
from recursive_horizons.nsc_covariant_identities import exact_checks
from recursive_horizons.nsc_covariant_measure import finite_cutoff_cylinder_stress
from recursive_horizons.nsc_finite_terms import response_matrix,response_probes
from recursive_horizons.nsc_shape_response import StaticAxialMetric, smooth_metric as staggered_metric, spin_shape_difference
from recursive_horizons.nsc_energy_transfer import regional_operators, vacuum, expectation
OUTPUT=ROOT/'results/nsc-9-covariant-source.json'


def native(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {key:native(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [native(item) for item in value]
    return value


def compare_record(expected,actual,path='$'):
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or expected.keys()!=actual.keys():raise RuntimeError(f'keys differ at {path}')
        for key in expected:compare_record(expected[key],actual[key],path+'/'+key)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(expected)!=len(actual):raise RuntimeError(f'list differs at {path}')
        for i,(a,b) in enumerate(zip(expected,actual)):compare_record(a,b,f'{path}/{i}')
    elif isinstance(expected,float):
        if isinstance(actual,bool) or not isinstance(actual,(float,int)) or not np.isclose(expected,actual,atol=2e-7,rtol=2e-8):
            raise RuntimeError(f'numeric field differs at {path}')
    elif type(expected) is not type(actual) or expected!=actual:raise RuntimeError(f'exact field differs at {path}')


def spectrum_controls():
    result=[]
    for eta,points in ((0.,17),(.5,16)):
        metric=smooth_metric(points,eta=eta,general=True)
        omega=.7;kappa=2
        plus=euclidean_operator(metric,omega,kappa)
        minus=euclidean_operator(metric,omega,-kappa)
        paired=np.block([[plus,np.zeros_like(plus)],[np.zeros_like(plus),minus]])
        gamma5=-np.kron(np.array([[0,1],[1,0]]),np.kron(SIGMA2,np.eye(points)))
        chirality=float(np.max(abs(paired@gamma5+gamma5@paired)))
        pq=radial_momentum(metric);f=np.diag(omega/metric.lapse);w=np.diag(kappa/metric.sphere_radius)
        square=(np.kron(np.eye(2),pq@pq+f@f+w@w)
                +1j*np.kron(SIGMA2,f@pq-pq@f)-1j*np.kron(SIGMA3,pq@w-w@pq))
        square_error=float(np.max(abs(plus@plus-square)))
        traces=[frequency_density(metric,o,k,2.,False)[0] for o,k in ((omega,kappa),(-omega,kappa),(omega,-kappa))]
        sigma=.15*np.cos(2*np.pi*metric.x/metric.length)+.07*np.sin(4*np.pi*metric.x/metric.length)
        s=np.diag(np.r_[sigma,sigma]);p=metric.momentum_matrix
        a=sigma/metric.radial_scale
        tangent=(-np.kron(SIGMA3,np.diag(omega*sigma/metric.lapse))
                 -.5*np.kron(SIGMA1,a[:,None]*p+p*a[None,:])
                 +np.kron(SIGMA2,np.diag(kappa*sigma/metric.sphere_radius)))
        sandwich=-.5*(s@plus+plus@s)
        assert chirality<1e-11 and square_error<1e-10 and np.ptp(traces)<1e-11
        result.append({'eta':eta,'points':points,'omega':omega,'kappa':kappa,
            'full_paired_chirality_residual':chirality,'finite_commutator_square_residual':square_error,
            'frequency_and_angular_sign_modulus_traces':traces,
            'finite_metric_vs_Weyl_sandwich_operator_difference_norm':float(np.linalg.norm(tangent-sandwich)),
            'physical_chirality_is_not_single_radial_grading':True})
    return result


def metric_variation_controls(metric,response):
    angle=2*np.pi*metric.x/metric.length;rows=[]
    for j,field in enumerate(('lapse','radial_scale','sphere_radius')):
        direction=getattr(metric,field)*(.3*np.cos(angle)+.2*np.sin(2*angle))
        analytic=float(np.dot(response['metric_gradients'][j],direction));controls=[]
        for step in (.002,.001):
            plus=replace(metric,**{field:getattr(metric,field)+step*direction})
            minus=replace(metric,**{field:getattr(metric,field)-step*direction})
            ep=cutoff_response(plus,angular_max=16,gradients=False)['energy']
            em=cutoff_response(minus,angular_max=16,gradients=False)['energy']
            fd=(ep-em)/(2*step)
            assert abs(fd-analytic)<1e-7
            controls.append({'step':step,'energy_plus':ep,'energy_minus':em,'finite_derivative':fd,'error':abs(fd-analytic)})
        rows.append({'field':field,'HF_derivative':analytic,'controls':controls})
    return rows


def schur_sign_controls():
    """A reciprocal, a frame change and a state current are different operations."""
    hp=np.diag([.8,-.8]);hc=1.2*SIGMA3+.25*SIGMA1
    b=.25*np.eye(2)+.12*SIGMA1+.04*SIGMA2
    z=.3+.4j;kc=z*np.eye(2)-hc
    sigma=b@np.linalg.solve(kc,b.conj().T)
    u=SIGMA1;changed_b=b@u.conj().T
    transported=changed_b@np.linalg.solve(u@kc@u.conj().T,changed_b.conj().T)
    frame_error=float(np.max(abs(transported-sigma)))
    joined=np.block([[hp,b],[b.conj().T,hc]])
    _,_,c=vacuum(joined)
    ops=regional_operators(joined,[1,1,0,0])
    current=expectation(c,ops['current_a'])
    assert frame_error<1e-12 and abs(current)<1e-12
    rows=[]
    for energy in (-2.,2.):
        child=energy*np.eye(2)-hc
        correction=-b@np.linalg.solve(child,b.conj().T)
        rows.append({'energy':energy,'child_inverse_eigenvalues':np.linalg.eigvalsh(np.linalg.inv(child)),
                     'Schur_correction_eigenvalues':np.linalg.eigvalsh(correction)})
    return {'noncommuting_complex_energy':{'real':z.real,'imag':z.imag},
            'self_energy_norm':float(np.linalg.norm(sigma)),
            'consistent_child_frame_change_residual':frame_error,
            'vacuum_current':current,'current_operator_norm':float(np.linalg.norm(ops['current_a'])),
            'real_energy_sign_controls':rows,
            'interpretation':'inverse is reciprocal, not sign reversal; Schur subtraction is elimination, not a state energy drain; physical expansion requires metric stress and dynamics',
            'cosmological_acceleration_inferred_from_Schur_sign':False}


def calculate():
    exact=exact_checks()
    cylinders=[]
    for eta,points in ((0.,25),(.5,24)):
        metric=cylinder_metric(points,eta=eta)
        r=cutoff_response(metric,angular_max=16)
        zero=cylinder_zero_winding_energy(angular_max=16)
        independent=finite_cutoff_cylinder_stress(4.,cutoff=2.,eta=eta)
        errors=[abs(r['energy']-ultrastatic_reference(metric,angular_max=16)),
                abs(r['energy']-zero-independent['energy']),
                max(abs(r['radial_null']-independent['axial_null_stress_integral']))]
        assert max(errors)<2e-10
        cylinders.append({'eta':eta,'points':points,'total_energy':r['energy'],
            'zero_winding_energy':zero,'compactification_energy':r['energy']-zero,
            'independent_compactification_energy':independent['energy'],
            'radial_null_stress':float(r['radial_null'][0]),
            'independent_radial_null_stress':independent['axial_null_stress_integral'],
            'ultrastatic_winding_and_null_errors':errors})

    spatial=[];responses={}
    for points in (16,24,32,48):
        metric=smooth_metric(points,general=True);r=cutoff_response(metric,angular_max=16)
        responses[points]=r
        spatial.append({'points':points,'energy':r['energy'],
            'neck_null_stress':float(r['radial_null'][points//2]),'wards':ward_summary(metric,r),
            'frequency_extent':r['frequency_extent'],
            'finite_angular_frequency_energy_tail_bound':r['finite_angular_frequency_energy_tail_bound']})
    assert spatial[-1]['wards']['maximum_radial_conservation_residual']<1e-9
    assert abs(spatial[-1]['energy']-spatial[-2]['energy'])<1e-9
    metric=smooth_metric(24,general=True);reference=responses[24]
    angular=[{'angular_max':k,**cutoff_response(metric,angular_max=k,gradients=False)} for k in (4,8,12,16,20)]
    assert abs(angular[-1]['energy']-angular[-2]['energy'])<2e-10
    frequency=[{'points':n,'energy':cutoff_response(metric,angular_max=16,frequency_points=n,gradients=False)['energy']} for n in (32,48,64)]
    assert abs(frequency[-1]['energy']-frequency[-2]['energy'])<1e-9
    extended=cutoff_response(metric,angular_max=16,frequency_points=64,extent_factor=1.25)
    assert abs(extended['energy']-reference['energy'])<1e-9
    sigma=.15*np.cos(2*np.pi*metric.x/metric.length)+.07*np.sin(4*np.pi*metric.x/metric.length)
    metric_weyl=float(np.sum(reference['metric_gradients']*np.array([metric.lapse,metric.radial_scale,metric.sphere_radius])*sigma))
    direct_weyl=cutoff_weyl_trace(metric,sigma)
    assert abs(metric_weyl-direct_weyl)<1e-9

    variation=metric_variation_controls(metric,reference)
    small=smooth_metric(16,general=True);angle=2*np.pi*small.x/small.length
    direction=np.array([.1*np.cos(angle),.07*np.sin(angle),.12*np.cos(2*angle)])
    differential=directional_variation(small,direction,angular_max=8)
    base=cutoff_response(small,angular_max=8,gradients=False)['energy'];hessian=[]
    for step in (.004,.002,.001):
        def changed(sign):return replace(small,**{field:getattr(small,field)*np.exp(sign*step*d)
                    for field,d in zip(('lapse','radial_scale','sphere_radius'),direction)})
        ep=cutoff_response(changed(1),angular_max=8,gradients=False)['energy']
        em=cutoff_response(changed(-1),angular_max=8,gradients=False)['energy']
        first=(ep-em)/(2*step);second=(ep+em-2*base)/step**2
        assert abs(first-differential['first'])<1e-7 and abs(second-differential['second'])<2e-7
        hessian.append({'step':step,'energy_plus':ep,'energy_minus':em,'finite_first':first,'finite_second':second})

    clock_metric=smooth_metric(32);clock_base=cutoff_response(clock_metric,angular_max=16)
    clock_changed=replace(clock_metric,lapse=1.2*clock_metric.lapse)
    clock_value=cutoff_response(clock_changed,angular_max=16)['energy']
    assert abs(clock_value-1.2*clock_base['energy'])<1e-9
    clock={'factor':1.2,'base_covariant_energy':clock_base['energy'],'changed_covariant_energy':clock_value,
           'changed_spatial_cutoff_proxy':spatial_cutoff_proxy(clock_changed),
           'interpretation':'fixed spatial-H cutoff is a different regulator; uniform clock scaling must multiply coordinate-time energy'}
    scaled=replace(metric,length=1.7*metric.length,sphere_radius=1.7*metric.sphere_radius)
    units=cutoff_response(scaled,cutoff=2/1.7,angular_max=16)
    assert abs(1.7*units['energy']-reference['energy'])<1e-9
    assert max(abs(units['radial_null']*1.7**4-reference['radial_null']))<1e-9

    heat_metric=smooth_metric(48,general=True);coefficients=heat_coefficients(heat_metric)
    heat=[]
    for time,angular_max in ((.2,24),(.1,32),(.05,48)):
        observed=heat_trace(heat_metric,time,angular_max=angular_max)
        predicted=(coefficients['a0']+time*coefficients['a2']+time*time*coefficients['a4'])/(4*np.pi*time)**2
        heat.append({'proper_time':time,'angular_max':angular_max,'trace_per_coordinate_time':observed,
                     'through_a4':predicted,'remainder':observed-predicted,'remainder_divided_by_time':(observed-predicted)/time})
    assert abs(heat[-1]['remainder'])<abs(heat[0]['remainder'])

    subtraction=[];m=smooth_metric(48);anomaly=local_trace_coefficient(m)
    for cutoff,angular_max in ((1.,12),(2.,20),(3.,28),(4.,36),(6.,48)):
        r=subtracted_response(m,cutoff,angular_max=angular_max,frequency_points=64)
        assert r['radial_null'][24]<0
        subtraction.append({'angular_max':angular_max,**r,
            'expected_subtraction_trace_limit':anomaly,
            'maximum_trace_limit_error':float(max(abs(r['trace']-anomaly)))})
    assert subtraction[-1]['maximum_trace_limit_error']<subtraction[1]['maximum_trace_limit_error']
    old=ultraviolet_subtraction(m,4.,1.);shifted=ultraviolet_subtraction(m,4.,np.e)
    shift=old['energy']-shifted['energy']
    assert abs(shift-old['remainder_normalization_log_derivative'])<1e-10
    existing=StaticAxialMetric(m.length,m.lapse,m.radial_scale,m.sphere_radius)
    probes=response_probes(existing);r=subtraction[3]
    source_projections=[float(np.sum(delta*r['metric_gradients'])) for _,delta in probes]

    # Independent prior staggered vacuum response, at a common geometric neck.
    pair=[]
    for eta,points in ((0.,49),(.5,48)):
        pmetric=smooth_metric(points,eta=eta)
        pr=cutoff_response(pmetric,cutoff=3.,angular_max=28,frequency_points=64)
        pair.append({'eta':eta,'points':points,'energy':pr['energy'],
                     'neck_null_stress':float(pr['radial_null'][points//2])})
    old_pair=[]
    for points in (128,256,512):
        previous=spin_shape_difference(staggered_metric(2.,points),8)
        old_pair.append({'points':points,'energy_difference':previous['energy_difference'],
                         'neck_null_difference':float(previous['axial_null_difference'][points//2])})
    extrapolated={key:(4*old_pair[-1][key]-old_pair[-2][key])/3 for key in ('energy_difference','neck_null_difference')}
    new_difference={'energy_difference':pair[0]['energy']-pair[1]['energy'],
                    'neck_null_difference':pair[0]['neck_null_stress']-pair[1]['neck_null_stress']}
    assert max(abs(new_difference[key]-extrapolated[key]) for key in new_difference)<2e-8

    source_paths=('src/recursive_horizons/nsc_covariant_operator.py','src/recursive_horizons/nsc_covariant_identities.py',
        'scripts/check_nsc_covariant_source.py','tests/test_nsc_covariant_source.py','tests/test_nsc_covariant_identities.py',
        'docs/nsc-covariant-source.md','src/recursive_horizons/nsc_regulated.py','src/recursive_horizons/nsc_finite_terms.py',
        'src/recursive_horizons/nsc_covariant_measure.py','src/recursive_horizons/nsc_shape_response.py',
        'src/recursive_horizons/nsc_energy_transfer.py')
    return native({'schema':'nsc-covariant-source-v1','artifact_id':'NSC-9-COVARIANT-SOURCE',
        'classification':'computed_static_four_dimensional_covariant_cutoff_modulus_and_explicit_UV_subtracted_metric_source',
        'source_hashes':{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in source_paths},
        'input_hashes':{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in (
            'results/nsc-4-covariant-measure.json','results/nsc-8-finite-terms.json')},
        'conventions':{'metric':'positive Euclidean N²dtau²+q²dx²+r²dOmega²; static continuation to signature+---',
            'measure':'chi=r sqrt(Nq) psi, flat dx dtau after angular separation',
            'operator':'D_omega=beta omega/N+i beta H0(q,r), full paired angular signs retained in counting',
            'numerical_reduction':'one two-component angular sign,4*kappa multiplicity; positive omega integral /pi; verified frequency/angular sign trace equality',
            'functional':'E_Lambda=.5 sum_kappa 2*kappa int_R domega/(2pi) Tr_(x,4) E1(D_omega²/Lambda²)',
            'state':'static vacuum on the declared compact axial P/AP domain; no trapped horizon continuation',
            'subtraction':'E_sub=E_Lambda-[a0 Lambda4/2+a2 Lambda2+a4 log(Lambda²/M²)]/(32pi²)',
            'scales':'Lambda is proper-time cutoff; M labels subtraction, not a fermion mass; neither selected physically',
            'phase':'modulus on invertible paired fibers; global phase/zero-mode prescription and common compensator remain open',
            'stress_trace':'rho-p_x-2p_perp=+Weyl variation density; opposite sign to conventions with deltaW=-int sigma Ttrace',
            'comparison':'all fields, float atol2e-7 rtol2e-8; keys integers booleans strings hashes exact; primitive control gates tighter',
            'finite_rank_normalization_not_integrated':True},
        'exact_four_dimensional_identities':exact,'spectrum_and_domain_controls':spectrum_controls(),
        'cylinder_controls':cylinders,'spatial_refinement':spatial,'angular_refinement':angular,
        'frequency_refinement':frequency,'extended_frequency_energy':extended['energy'],
        'local_Weyl_control':{'metric_variation':metric_weyl,'raw_cutoff_trace':direct_weyl},
        'metric_HF_controls':variation,'independent_divided_difference_variation':differential,
        'hessian_finite_differences':hessian,'hessian_base_energy':base,'uniform_clock_control':clock,
        'unit_scaling':{'factor':1.7,'scaled_energy':units['energy'],'expected_scaled_energy':reference['energy']/1.7},
        'full_nonconstant_metric_response':{'x':smooth_metric(48,general=True).x,**responses[48]},
        'ultraviolet_heat_coefficients':coefficients,'heat_asymptotic_controls':heat,
        'subtracted_determinant_sources':subtraction,
        'normalization_change':{'log_M_change':1.,'energy_shift':shift,'predicted_shift':old['remainder_normalization_log_derivative']},
        'Schur_sign_and_energy_controls':schur_sign_controls(),
        'metric_matching':{'cutoff':4.,'M':1.,'points':48,'probe_names':[name for name,_ in probes],
            'computed_determinant_source_projections':source_projections,
            'local_response_matrix':response_matrix(existing),
            'equation':'b_A_determinant + J_Ai c_i + b_A_compensator_link_state = 0',
            'coefficients_fitted_to_imposed_neck':False},
        'same_geometry_spin_structure_bridge':{'covariant_cutoff':3.,'covariant_angular_max':28,
            'covariant_states':pair,'covariant_difference':new_difference,
            'independent_staggered_vacuum_differences':old_pair,'second_order_extrapolation':extrapolated,
            'extrapolation_is_a_rigorous_error_bound':False},
        'nonclaims':{'complete_invariant_functional_fixed':False,'finite_coefficients_selected_by_anomaly':False,
            'subtracted_null_sign_is_scheme_independent':False,'self_sourced_stationary_metric':False,
            'retarded_metric_response_or_in_in_dynamics_computed':False,'physical_Phi_or_Omega_derived':False,
            'new_to_world_heat_coefficients':False,'full_physical_closure':False},'terminal':True})


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite record')
    record=calculate()
    if args.check:compare_record(json.loads(OUTPUT.read_text()),record);print('Covariant Dirac source, UV coefficients and all metric controls reproduced; every field checked.')
    elif args.output:
        with args.output.open('x') as stream:json.dump(record,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        print(args.output)
    else:print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
