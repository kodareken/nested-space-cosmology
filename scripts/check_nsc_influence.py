#!/usr/bin/env python3
"""Reproduce normalized Dirac influence, causal response and geometric noise."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_covariant_operator import smooth_metric,euclidean_operator,SIGMA3
from recursive_horizons.nsc_influence import (
    canonical_hamiltonian,geometric_vertex,ground_covariance,thermal_covariance,
    influence,mean,connected,noise_kernel,retarded_response,gaussian_smeared_vertex,
    response_spectrum,zero_frequency_noise,pulse_unitary,branch_unitary,
    mixed_response_from_influence,static_susceptibility,with_real_branch_action,
    fock_annihilators,second_quantize,gaussian_fock_state,
)
from recursive_horizons.nsc_vacuum_work import leading_response
OUTPUT=ROOT/'results/nsc-10-influence.json'


def native(value):
    if isinstance(value,np.ndarray):return native(value.tolist())
    if isinstance(value,np.generic):return native(value.item())
    if isinstance(value,complex):return {'real':float(value.real),'imag':float(value.imag)}
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
        if isinstance(actual,bool) or not isinstance(actual,(float,int)) or not np.isclose(expected,actual,rtol=2e-8,atol=2e-7):
            raise RuntimeError(f'number differs at {path}')
    elif type(expected) is not type(actual) or expected!=actual:raise RuntimeError(f'exact field differs at {path}')


def fock_controls():
    rng=np.random.default_rng(419);a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));h=(a+a.conj().T)/3
    a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));v=(a+a.conj().T)/3
    _,u=np.linalg.eigh(v);c=(u*np.array([.2,.6,.8])[None,:])@u.conj().T
    ann=fock_annihilators(3);rho=gaussian_fock_state(c,ann)
    restored=np.array([[np.trace(rho@ann[j].conj().T@ann[i]) for j in range(3)] for i in range(3)])
    covariance_error=float(np.max(abs(restored-c)))
    ph=expm(-.4j*(h+.23*v));mh=expm(-.4j*(h-.17*v))
    value=influence(c,ph,mh)
    direct=np.trace(expm(-.4j*second_quantize(h+.23*v,ann))@rho@expm(.4j*second_quantize(h-.17*v,ann)))
    error=float(abs(value['amplitude']-direct))
    # Pure-state limit: an occupied orbital is constructed by a Fock creation operator.
    orbital=u[:,1];cp=np.outer(orbital,orbital.conj())
    creation=sum(orbital[i]*ann[i].conj().T for i in range(3))
    state=creation@np.eye(8)[:,0];pure_rho=np.outer(state,state.conj())
    pure_direct=np.trace(expm(-.4j*second_quantize(h+.23*v,ann))@pure_rho@expm(.4j*second_quantize(h-.17*v,ann)))
    pure_det=influence(cp,ph,mh)['amplitude']
    assert max(covariance_error,error,abs(pure_direct-pure_det))<1e-12
    return {'one_particle_modes':3,'Fock_dimension':8,'occupations':[.2,.6,.8],
            'covariance_reconstruction_error':covariance_error,'determinant':value['amplitude'],
            'direct_Fock_trace':direct,'mixed_state_error':error,
            'pure_orbital_determinant':pure_det,'pure_orbital_Fock_trace':pure_direct}


def calculate():
    fock=fock_controls();metric=smooth_metric(16);h=canonical_hamiltonian(metric);v=geometric_vertex(metric)
    e,u,c=ground_covariance(h);identity=np.eye(len(h));free=expm(-.3j*h)
    same=influence(c,free,free)
    assert abs(same['amplitude']-1)<1e-12
    general=smooth_metric(16,general=True);hg=canonical_hamiltonian(general)
    ni=np.diag(np.r_[general.lapse,general.lapse]**-.5);beta=np.kron(SIGMA3,np.eye(general.points))
    covariant_error=float(np.max(abs(euclidean_operator(general,.7)-ni@beta@(.7*np.eye(len(hg))+1j*hg)@ni)))
    assert covariant_error<1e-11

    smeared=gaussian_smeared_vertex(h,v,.5);mu=mean(c,smeared).real;variance=connected(c,smeared,smeared).real
    derivatives=[]
    for delta in (.008,.004,.002):
        plus=influence(c,expm(-.5j*delta*smeared),expm(.5j*delta*smeared))
        minus=influence(c,expm(.5j*delta*smeared),expm(-.5j*delta*smeared))
        first=(plus['principal_action']-minus['principal_action'])/(2*delta)
        second=(plus['principal_action']+minus['principal_action'])/delta**2
        assert abs(first.real+mu)<3e-7 and abs(second.imag-variance)<1e-7
        assert abs(plus['amplitude'].conjugate()-minus['amplitude'])<1e-12
        derivatives.append({'delta':delta,'first_difference_derivative':first,
                            'second_difference_derivative':second,'overlap_modulus':abs(plus['amplitude'])})
    times=[0.,.2,.6,1.];noise=noise_kernel(h,c,v,times);eigen=np.linalg.eigvalsh(noise)
    assert min(eigen)>-1e-10
    causal=[]
    for source in (.2,1.1):
        expected=-retarded_response(h,c,v,v,.7,source)
        observed=mixed_response_from_influence(h,c,v,v,.7,source)
        assert abs(observed-expected)<2e-7
        causal.append({'observation_time':.7,'source_time':source,'mixed_CTP_derivative':observed,'minus_retarded_response':expected})
    susceptibility=static_susceptibility(h,v);ground=float(sum(e[e<0]));static=[]
    for step in (.004,.002):
        ep=np.linalg.eigvalsh(h+step*v);em=np.linalg.eigvalsh(h-step*v)
        second=(sum(ep[ep<0])+sum(em[em<0])-2*ground)/step**2
        assert abs(second-susceptibility)<2e-7
        static.append({'step':step,'ground_energy_second_difference':float(second)})

    thermal=thermal_covariance(h,1.7);spectral=response_spectrum(h,thermal,v)
    residual=max(abs(row['noise_plus']-row['minus_im_chi_plus']/np.tanh(.85*row['frequency'])) for row in spectral)
    assert residual<1e-12
    resolution=[]
    for points in (12,16,24,32):
        m=smooth_metric(points);hh=canonical_hamiltonian(m);vv=geometric_vertex(m);_,_,cc=ground_covariance(hh)
        sm=gaussian_smeared_vertex(hh,vv,.5);prediction=leading_response(hh,vv)
        var=connected(cc,sm,sm).real
        assert abs(var-prediction['pair_number_coefficient'])<1e-12
        resolution.append({'points':points,'bare_mean_vertex':mean(cc,vv).real,
            'equal_time_variance':connected(cc,vv,vv).real,'smeared_variance':var,
            'pair_number_coefficient':prediction['pair_number_coefficient'],
            'energy_coefficient':prediction['energy_coefficient']})
    assert abs(resolution[-1]['smeared_variance']-resolution[-2]['smeared_variance'])<1e-10

    histories=[]
    for steps in (64,128,256):
        plus=pulse_unitary(h,v,.02,steps=steps)
        zero=np.eye(len(h));value=influence(c,plus,zero);evolved=plus@c@plus.conj().T
        particles=float(np.trace((identity-c)@evolved).real)
        holes=float(np.trace(c@(identity-evolved)).real)
        energy=float(np.trace(h@(evolved-c)).real)
        occ=u[:,e<0];overlap=np.linalg.det(occ.conj().T@plus@occ)
        assert abs(overlap-value['amplitude'])<1e-11 and abs(particles-holes)<1e-10
        assert abs(value['amplitude'])<=1+1e-10 and energy>0
        histories.append({'steps':steps,'epsilon':.02,'width':.5,'influence':value,
            'Slater_overlap':overlap,'created_particles':particles,'created_holes':holes,
            'excitation_energy':energy,'unitarity_error':float(np.max(abs(plus.conj().T@plus-identity)))})
    # The exact pulse is time ordered. Its leading imaginary curvature is the
    # smeared covariance, but exp(-i epsilon V_f) is not the complete pulse U.
    small_histories=[]
    for epsilon in (.02,.01):
        plus=pulse_unitary(h,v,epsilon,steps=256);minus=pulse_unitary(h,v,-epsilon,steps=256)
        value=influence(c,plus,minus)
        coefficient=-value['log_modulus']/(2*epsilon**2)
        small_histories.append({'branch_amplitudes':[epsilon,-epsilon],
                                'negative_log_modulus_over_2epsilon_squared':coefficient})
    assert abs(small_histories[-1]['negative_log_modulus_over_2epsilon_squared']-variance)<1e-5

    phase=.37
    original=histories[-1]['influence']['amplitude']
    shifted=with_real_branch_action(original,phase*.02,0.)
    same_shift=with_real_branch_action(1.,phase*.02,phase*.02)
    assert abs(abs(shifted)-abs(original))<1e-14 and same_shift==1
    zero=influence(np.array([[.5]]),np.array([[-1.]]),np.array([[1.]]))
    assert zero['amplitude']==0 and zero['principal_action'] is None
    return native({'schema':'nsc-influence-v1','artifact_id':'NSC-10-INFLUENCE',
        'classification':'normalized_finite_Gaussian_Dirac_state_functional_with_causal_response_and_smeared_geometric_noise',
        'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'src/recursive_horizons/nsc_influence.py','scripts/check_nsc_influence.py','tests/test_nsc_influence.py',
            'docs/nsc-influence.md','src/recursive_horizons/nsc_covariant_operator.py',
            'src/recursive_horizons/nsc_vacuum_work.py')},
        'input_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'results/nsc-9-covariant-source.json','results/nsc-6-vacuum-work.json')},
        'conventions':{'units':'hbar=c=a=1; one representative kappa=1 channel unless explicitly stated',
            'geometry':'smooth R=2 AP axial cell, canonical H from the NSC9 covariant operator',
            'covariance':'Cij=<c_j dagger c_i>, 0<=C<=I',
            'Hamiltonian_source':'H_int=+J(t)V, radius r=r0/(1+J s), s from the same smooth profile',
            'influence':'Z=det(I-C+C U_minus dagger U_plus), Gamma_IF=-i logZ',
            'first_difference_derivative':'-<V>',
            'retarded':'chi_AB=-i theta(t-s) Tr C[A(t),B(s)], Fourier exp(+i omega t)',
            'noise':'Re Tr C A(t)(I-C)B(s)',
            'positive_frequency_FDT':'N(omega)=-coth(beta omega/2) Im chi_R(omega); elastic zero-frequency part separately retained',
            'branch':'principal local log branch near equal histories; zero overlap valid but log undefined',
            'cutoff':'finite spatial regulator; bare mean/equal-time noise are not renormalized sources',
            'pulse':'interaction-picture unitary for stationary initial vacuum, window +/-7 width; smooth-history time refinement separate',
            'probe':'source impulses test functional derivatives; not finite delta-function geometries',
            'comparison':'all fields; float atol2e-7 rtol2e-8; exact remaining types/keys/hashes'},
        'Fock_trace_control':fock,'equal_history':same,'covariant_to_canonical_operator_residual':covariant_error,
        'smeared_generator':{'width':.5,'mean':mu,'variance':variance,'local_log_derivatives':derivatives},
        'noise_matrix':{'times':times,'values':noise,'eigenvalues':eigen},'causal_contour_controls':causal,
        'static_susceptibility':{'retarded_zero_frequency':susceptibility,'finite_energy_variations':static},
        'thermal_FDT':{'inverse_temperature':1.7,'maximum_cross_multiplied_residual':residual,
            'zero_frequency_elastic_noise':zero_frequency_noise(h,thermal,v),'positive_frequency_weights':spectral},
        'spatial_and_time_smearing_refinement':resolution,'time_ordered_pulse_histories':histories,
        'opposite_pulse_small_amplitude_controls':small_histories,
        'local_real_phase_control':{'alpha':phase,'original_overlap':original,'shifted_overlap':shifted,
            'equal_history_after_phase':same_shift,'mean_source_shift':'alpha for F(J)=alpha J; illustrative, not selected physical coefficient'},
        'zero_overlap_control':{'amplitude':zero['amplitude'],'log_action_defined':False},
        'nonclaims':{'full_covariant_renormalized_in_in_action':False,'metric_constraints_or_graviton_health':False,
            'absolute_source_fixed_by_Z_equal_one':False,'Born_rule_derived':False,'physical_self_sourcing':False,
            'invariant_finite_completion_fixed':False,'new_to_world_Gaussian_trace_identity':False},'terminal':True})


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite record')
    record=calculate()
    if args.check:compare_record(json.loads(OUTPUT.read_text()),record);print('Normalized Dirac influence, causal response and geometric noise reproduced; every field checked.')
    elif args.output:
        with args.output.open('x') as stream:json.dump(record,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        print(args.output)
    else:print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
