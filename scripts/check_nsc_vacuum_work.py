#!/usr/bin/env python3
"""Verify geometric vacuum excitation, clock control and supplied-work balance."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_shape_response import smooth_metric, StaticAxialMetric
from recursive_horizons.nsc_energy_transfer import dirac_hamiltonian
from recursive_horizons.nsc_vacuum_work import radius_vertex, pulse_profile, leading_response, evolve_vacuum_pulse
OUTPUT = ROOT/'results/nsc-6-vacuum-work.json'


def compare_record(expected, observed, path=''):
    if isinstance(expected,bool) or isinstance(observed,bool):
        if expected is not observed: raise RuntimeError(f'boolean mismatch at {path}')
    elif isinstance(expected,int) and isinstance(observed,int):
        if expected != observed: raise RuntimeError(f'integer mismatch at {path}')
    elif isinstance(expected,(int,float)) and isinstance(observed,(int,float)):
        if not np.isclose(expected,observed,rtol=1e-7,atol=5e-10):
            raise RuntimeError(f'numeric mismatch at {path}')
    elif isinstance(expected,dict) and isinstance(observed,dict):
        if expected.keys()!=observed.keys(): raise RuntimeError(f'keys differ at {path}')
        for k in expected: compare_record(expected[k],observed[k],path+'/'+k)
    elif isinstance(expected,list) and isinstance(observed,list):
        if len(expected)!=len(observed): raise RuntimeError(f'length mismatch at {path}')
        for i,(a,b) in enumerate(zip(expected,observed)): compare_record(a,b,f'{path}/{i}')
    elif expected != observed: raise RuntimeError(f'value mismatch at {path}')


def calculate():
    t,x = sp.symbols('t x',real=True)
    r,u = sp.Function('r')(t,x),sp.Function('u')(t,x)
    exact = {str(variable): str(sp.simplify(r*(sp.diff(u/r,variable)+sp.diff(r,variable)/r*(u/r))-sp.diff(u,variable)))
             for variable in (t,x)}
    assert set(exact.values()) == {'0'}
    grids=[]
    for points in (32,64,128,256):
        m=smooth_metric(2.,points)
        h,v=dirac_hamiltonian(m),radius_vertex(m)
        response=leading_response(h,v)
        eps=.13
        changed=replace(m,sphere_radius=m.sphere_radius/(1+eps*pulse_profile(m)))
        affine_error=float(np.max(abs(dirac_hamiltonian(changed)-h-eps*v)))
        assert affine_error<1e-12
        full=sum(4*k*leading_response(dirac_hamiltonian(m,k),radius_vertex(m,k))['energy_coefficient'] for k in range(1,9))
        grids.append({'points':points,'response':response,'full_angular_through_kappa8_energy_coefficient':full,
                      'exact_geometric_affinity_residual':affine_error})

    m=smooth_metric(2.,128)
    angular=[]
    cumulative=0.
    for k in range(1,9):
        response=leading_response(dirac_hamiltonian(m,k),radius_vertex(m,k))
        cumulative+=4*k*response['energy_coefficient']
        angular.append({'kappa':k,'multiplicity':4*k,'one_channel':response,'cumulative_energy_coefficient':cumulative})
    assert angular[-1]['cumulative_energy_coefficient']-angular[-3]['cumulative_energy_coefficient']<1e-8
    widths=[{'width':width,**leading_response(dirac_hamiltonian(m),radius_vertex(m),width)} for width in (.25,.5,1.,2.)]
    clock=leading_response(dirac_hamiltonian(m),dirac_hamiltonian(m))
    assert clock['energy_coefficient']<1e-25

    # Fourier-symbol control on an independent constant-radius geometry.
    uniform=[]
    for eta in (0.,.5):
        cylinder=StaticAxialMetric(4.,np.ones(64),np.ones(64),np.ones(64))
        h=dirac_hamiltonian(cylinder,1,eta)
        vertex=radius_vertex(cylinder,1,eta,np.ones(64))
        matrix=leading_response(h,vertex)
        phase=2*np.pi*(np.arange(64)+eta)/64
        p=2*np.sin(phase/2)/cylinder.spacing
        mass=np.cos(phase/2)
        e=np.sqrt(p*p+mass*mass)
        transition=mass*mass*p*p/(e*e)
        weight=transition*(2*np.pi*.5**2)*np.exp(-e*e)
        residual=abs(matrix['energy_coefficient']-float(np.sum(2*e*weight)))
        assert residual<1e-12
        uniform.append({'eta':eta,'matrix':matrix,'independent_Fourier_energy_coefficient':float(np.sum(2*e*weight)),
                        'energy_residual':residual})

    evolved=[]
    for epsilon in (.04,.02,.01,-.04,-.02,-.01):
        row=evolve_vacuum_pulse(smooth_metric(2.,32),epsilon)
        assert row['work_balance_residual']<5e-10
        assert row['regional_balance_residual']<5e-10
        assert row['orthonormality_residual']<3e-12
        assert abs(row['pair_number']-row['hole_number'])<3e-11
        evolved.append(row)
    coefficient=grids[0]['response']['energy_coefficient']
    symmetric=[]
    for i in range(3):
        epsilon=evolved[i]['epsilon']
        value=(evolved[i]['external_work']+evolved[i+3]['external_work'])/(2*epsilon**2)
        symmetric.append({'epsilon':epsilon,'even_work_coefficient':value,'difference_from_leading':abs(value-coefficient)})
    assert symmetric[-1]['difference_from_leading']<symmetric[0]['difference_from_leading']/10
    fine=evolve_vacuum_pulse(smooth_metric(2.,32),.02,rtol=2e-13)
    assert abs(fine['external_work']-evolved[1]['external_work'])<1e-10
    assert fine['work_balance_residual']<1e-11
    assert fine['regional_balance_residual']<1e-11
    return {'schema':'nsc-vacuum-work-v1','artifact_id':'NSC-6-VACUUM-WORK',
        'classification':'finite_and_perturbative_Dirac_vacuum_response_to_prescribed_geometric_work_not_self_sourced_cosmology',
        'source_hashes':{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in (
            'src/recursive_horizons/nsc_vacuum_work.py','src/recursive_horizons/nsc_energy_transfer.py',
            'src/recursive_horizons/nsc_shape_response.py','scripts/check_nsc_vacuum_work.py')},
        'conventions':{'metric':'ds²=dt²-dx²-r(t,x)² dOmega², compact axial AP circle, original a=1,R=2',
            'radius':'r(t,x)=r0(x)/(1+epsilon*s(x)*g(t)); no added scalar matter',
            'spatial_profile':'s=exp(2*(cos(2*pi*(x+0.6)/L)-1))',
            'pulse':'g(t)=exp(-t²/(2*width²)), prescribed rather than action-derived',
            'state':'initial static vacuum; final particles defined relative to the same static H0',
            'time_reduction':'u=r(t,x)*psi removes both temporal and spatial radius connections',
            'units':'hbar=c=a_throat=1, width in a_throat, energy in 1/a_throat',
            'full_angular_multiplier':'4*kappa per radial representative; energy already counts particle and hole',
            'perturbative_domain':'coefficients multiply epsilon², generic omitted order epsilon³; sign average removes odd powers',
            'real_time_domain':'one kappa=1 representative, 32 axial points, interaction-picture DOP853, endpoints +/-7*width',
            'regional_ledger':'x>=0 half-link energy; integrated inflow includes throat and compact seam together; explicit geometric work is separate',
            'comparison':'all keys and integers exact, floats rtol=1e-7,atol=5e-10'},
        'exact_reduction_residuals':exact,
        'leading_identity':'P_nm/epsilon²=2*pi*width²*|V_nm|²*exp[-(E_n-E_m)²*width²]; DeltaE=sum_(n>0,m<0)(E_n-E_m)P_nm>=0',
        'work_identity':'d< H(t) >/dt=Tr(C dot H); for equal endpoint Hamiltonians, supplied work becomes excitation energy',
        'spatial_refinement':grids,'angular_refinement':angular,'duration_controls':widths,
        'pure_clock_control':clock,'independent_constant_geometry_controls':uniform,
        'finite_real_time_controls':evolved,'even_perturbation_controls':symmetric,'tighter_integrator_control':fine,
        'nonclaims':{'external_geometry_work_source_solved':False,'sustained_cosmological_Q_derived':False,
            'absolute_vacuum_stress_derived':False,'baryonic_dust_produced':False,
            'stationary_recursive_solution':False,'infinite_angular_remainder_bound_proved':False,
            'new_to_world_particle_creation_mechanism':False},'terminal':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists(): raise FileExistsError('refusing to overwrite result')
    record=calculate()
    if args.check:
        compare_record(json.loads(OUTPUT.read_text()),record)
        print('Geometric vacuum excitation and external-work balance reproduced; every field checked.')
    elif args.output:
        with args.output.open('x') as stream:
            json.dump(record,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        print(args.output)
    else: print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__': main()
