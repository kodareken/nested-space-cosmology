#!/usr/bin/env python3
"""Reproduce the smooth-geometry Euclidean/retarded radius-response match."""
from functools import lru_cache
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_covariant_operator import smooth_metric, directional_variation
from recursive_horizons.nsc_response_matching import (
    canonical_data, radius_path, radius_profile, static_proper_time,
    retarded_susceptibility, euclidean_kernel, rational_frequency_bubble, abel_extrapolation)
OUTPUT=ROOT/'results/nsc-11-response-matching.json'


def native(v):
    if isinstance(v,dict):return {k:native(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [native(x) for x in v]
    if isinstance(v,np.ndarray):return native(v.tolist())
    if isinstance(v,np.generic):return v.item()
    return v


def compare_record(expected,actual,path='$'):
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or expected.keys()!=actual.keys():raise RuntimeError('keys differ at '+path)
        for k in expected:compare_record(expected[k],actual[k],path+'/'+k)
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(expected)!=len(actual):raise RuntimeError('list differs at '+path)
        for i,(x,y) in enumerate(zip(expected,actual)):compare_record(x,y,f'{path}/{i}')
    elif isinstance(expected,float):
        if isinstance(actual,bool) or not isinstance(actual,(float,int)) or not np.isclose(expected,actual,atol=2e-8,rtol=2e-8):
            raise RuntimeError(f'numeric field differs at {path}: {expected} != {actual}')
    elif type(expected) is not type(actual) or expected!=actual:raise RuntimeError('exact field differs at '+path)


@lru_cache(None)
def metric(n,eta):return smooth_metric(n,eta=eta)


@lru_cache(None)
def data(n,eta,kappa):return canonical_data(metric(n,eta),kappa)


@lru_cache(None)
def euclidean(n,eta,kappa,nu,cutoff,nf):
    return euclidean_kernel(metric(n,eta),nu,cutoff,kappa,nf)


def relative(n,kappa,nu,cutoff,nf):
    return 4*kappa*(euclidean(n+1,0.,kappa,nu,cutoff,nf)-euclidean(n,.5,kappa,nu,cutoff,nf))


def calculate():
    small=metric(12,.5);d=data(12,.5,1)
    finite=[]
    for nu in (0.,.6,1.4):
        direct=rational_frequency_bubble(small,nu)
        ret=retarded_susceptibility(d,1j*nu).real
        finite.append({'nu':nu,'direct_infinite_frequency_integral':direct,
                       'retarded_at_i_nu':ret,'difference':direct-ret})
    static=static_proper_time(d,2.)
    contact_controls=[]
    for coordinate,key in [('inverse_radius','inverse_radius_hessian'),('log_radius','log_radius_hessian')]:
        rows=[]
        for step in (.008,.004):
            ep=static_proper_time(canonical_data(radius_path(small,step,coordinate)),2.)['energy']
            em=static_proper_time(canonical_data(radius_path(small,-step,coordinate)),2.)['energy']
            fd=(ep-2*static['energy']+em)/step**2
            rows.append({'step':step,'energy_second_difference':fd,'minus_spectral_hessian':fd-static[key]})
        contact_controls.append({'coordinate':coordinate,'hessian':static[key],'differences':rows})
    log_direction=np.array([np.zeros(12),np.zeros(12),radius_profile(small)])
    old=directional_variation(small,log_direction,cutoff=2.,angular_max=1,frequency_points=64)
    old_hessian=old['second']
    covariant_match={'full_kappa1_weight':4,'frequency_operator_hessian':old_hessian,
                     'canonical_log_radius_hessian_times4':4*static['log_radius_hessian'],
                     'difference':old_hessian-4*static['log_radius_hessian']}
    spatial=[]
    for n in (16,24,32,48):
        spatial.append({'P_points':n+1,'AP_points':n,
                        'nu1p4_kappa1_difference':relative(n,1,1.4,4.,96)})
    frequency=[]
    for nf in (32,48,64,96):
        frequency.append({'nodes':nf,'nu1p4_kappa1_difference':relative(32,1,1.4,4.,nf)})
    cutoff=[]
    for scale in (1.,2.,3.,4.,6.):
        cutoff.append({'proper_time_cutoff':scale,'nu1p4_kappa1_difference':relative(48,1,1.4,scale,96)})
    angular=[]
    for nu in (0.,.6,1.4):
        total=0.;rows=[]
        for kappa in range(1,9):
            value=relative(32,kappa,nu,4.,96)
            total+=value
            rows.append({'kappa':kappa,'weighted_difference':value,'cumulative':total})
        angular.append({'nu':nu,'sectors':rows})
    abel=[]
    for n,scale in ((128,20.),(192,32.),(256,48.)):
        for nu in (0.,.6,1.4):
            rows=[];total=0.
            for kappa in range(1,9):
                row=abel_extrapolation(data(n+1,0.,kappa),data(n,.5,kappa),nu,scale)
                row['kappa']=kappa;row['weight']=4*kappa
                total+=4*kappa*row['extrapolated'];rows.append(row)
            target=next(x['sectors'][-1]['cumulative'] for x in angular if x['nu']==nu)
            abel.append({'P_points':n+1,'AP_points':n,'nu':nu,'sectors':rows,
                         'total':total,'minus_euclidean':total-target})
    # Spatial cutoff and Abel cutoff are refined separately for the dominant channel.
    abel_spatial=[]
    for n in (128,192,256):
        row=abel_extrapolation(data(n+1,0.,1),data(n,.5,1),1.4,32.)
        abel_spatial.append({'P_points':n+1,'AP_points':n,**row})
    contact_pair=[]
    for eta,n in ((0.,49),(.5,48)):
        s=static_proper_time(data(n,eta,1),4.)
        contact_pair.append({'spin_structure':'P' if eta==0 else 'AP',**s})
    relative_contact=contact_pair[0]['log_radius_contact']-contact_pair[1]['log_radius_contact']
    finite_cutoff={'proper_time_kappa1_static':static['inverse_radius_hessian'],
                   'unregulated_fixed_grid_vacuum_kubo':finite[0]['retarded_at_i_nu'],
                   'equal_at_finite_cutoff':False,
                   'meaning':'Proper-time regularization of the action is not the bare finite-grid Kubo kernel or an Abel transition cutoff.'}
    gates={
        'finite_frequency_loop_matches_retarded':max(abs(x['difference']) for x in finite)<1e-11,
        'same_covariant_static_hessian':abs(covariant_match['difference'])<2e-8,
        'independent_energy_hessians':max(abs(x['differences'][-1]['minus_spectral_hessian']) for x in contact_controls)<2e-6,
        'spatial_refinement':abs(spatial[-1]['nu1p4_kappa1_difference']-spatial[-2]['nu1p4_kappa1_difference'])<2e-8,
        'frequency_refinement':abs(frequency[-1]['nu1p4_kappa1_difference']-frequency[-2]['nu1p4_kappa1_difference'])<2e-8,
        'angular_last_increment_small':max(abs(x['sectors'][-1]['weighted_difference']) for x in angular)<1e-7,
        'relative_euclidean_and_Abel_limit_match':max(abs(x['minus_euclidean']) for x in abel if x['AP_points']==256)<5e-8,
        'state_dependent_coordinate_contact_retained':abs(relative_contact)>1e-4,
    }
    if not all(gates.values()):raise RuntimeError(f'matching gate failed: {gates}')
    sources=('scripts/check_nsc_response_matching.py','src/recursive_horizons/nsc_response_matching.py',
             'src/recursive_horizons/nsc_covariant_operator.py','src/recursive_horizons/nsc_covariant_identities.py',
             'src/recursive_horizons/nsc_influence.py','tests/test_nsc_response_matching.py','docs/nsc-response-matching.md')
    inputs=('results/nsc-9-covariant-source.json','results/nsc-10-influence.json','results/nsc-10-measure-normalization.json')
    return native({'schema':'nsc-response-matching-v1','artifact_id':'NSC-11-RESPONSE-MATCHING',
        'classification':'relative_Euclidean_and_retarded_radius_response_matching_on_smooth_ultrastatic_cell',
        'conventions':{'geometry':'R2 smooth circle, a=1,N=q=1; P and AP spin structures on symmetric Fourier grids',
            'state':'negative-energy Gaussian vacuum in each spin structure, not a dynamically selected state',
            'coordinate':'r=r0/(1+J s), s=exp[2(cos(2pi(x+0.6)/L)-1)]; affine Dirac perturbation',
            'units':'hbar=c=a=1; energies and frequencies in inverse a; J dimensionless',
            'angular_counting':'representative2N block times4kappa; finite angular truncation explicit',
            'retarded':'H_int=+J V, Fourier exp(+i z t), analytic for Im z>0',
            'relative':'P minus AP on identical local geometry; local state-independent finite metric terms cancel in continuum',
            'contact':'log r coordinate has <H_second>; it is state-dependent and need not cancel in P-AP',
            'cutoffs':'proper-time Lambda and positive Abel transition weight are distinct finite prescriptions; matching tests cutoff-removal limit'},
        'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
        'input_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        'finite_grid_frequency_identity':finite,'static_spectral':static,'independent_metric_differences':contact_controls,
        'covariant_static_owner_match':covariant_match,'finite_cutoff_distinction':finite_cutoff,
        'spatial_resolution':spatial,'frequency_resolution':frequency,'proper_time_cutoff_resolution':cutoff,
        'angular_resolution':angular,'Abel_resolution':abel,'Abel_fixed_cutoff_spatial_resolution':abel_spatial,
        'contact_spin_structures':contact_pair,'one_channel_relative_log_radius_contact':relative_contact,
        'extrapolation':'quadratic Richardson in1/K^2; cancels first two analytic Abel terms; measured convergence, no proved remainder bound',
        'comparison':'every field; floats atol2e-8 rtol2e-8, exact keys types hashes; internal identity gates tighter',
        'gate':gates,'nonclaims':{'complete_finite_cutoff_causal_functional':False,'local_finite_coefficients_fixed':False,
            'preferred_spin_structure_selected':False,'metric_hessian_constraints_or_health':False,
            'self_sourced_solution':False,'continuum_remainder_bound':False,'new_to_world_Kubo_identity':False},'terminal':True})


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite record')
    record=calculate()
    if args.check:compare_record(json.loads(OUTPUT.read_text()),record);print('Euclidean/retarded geometric response and coordinate contacts reproduced; all fields checked.')
    elif args.output:
        with args.output.open('x') as f:json.dump(record,f,indent=2,sort_keys=True,allow_nan=False);f.write('\n')
        print(args.output)
    else:print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
