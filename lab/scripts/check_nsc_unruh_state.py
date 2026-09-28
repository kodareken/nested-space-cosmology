#!/usr/bin/env python3
"""Reproduce only the parent-matched massless source and its interface checks."""
import argparse
import json
from math import pi
from pathlib import Path
import sys

import numpy as np
from scipy.special import zeta
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_unruh_state import (
    ParentDirac,horizon_covariance,horizon_frame,rindler_ratio_control,
    parent_angular_source,sixth_order_angular_coefficient,thermal_phase_omission_bound,
)
from recursive_horizons.nsc_angular_stress import physical_source
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT=ROOT/'results/development/unruh-state.json'
SOURCES=('scripts/check_nsc_unruh_state.py','src/recursive_horizons/nsc_unruh_state.py',
         'docs/nsc-unruh-state.md','src/recursive_horizons/nsc_gauge_source.py',
         'scripts/check_nsc_compact_boundary_action.py','scripts/check_nsc_compact_casimir.py',
         'scripts/check_nsc_vacuum_charge_matching.py')
INPUTS=('results/development/horizon-source.json','results/development/angular-stress.json',
        'results/development/gauge-source.json')
EDGES=(0.,.25,1.,4.,16.,64.,128.,256.)
POINTS=(.1,.2,.5,1.,pi/2)
COMPONENTS=('rho','p_parallel','p_sphere','radial_null','T_hatT_hatz')
STATE='Unruh_boundary_candidate'
SPECS={
    'n16':dict(angular_max=16,points_per_interval=32,frequency_edges=EDGES,steps=24000),
    'n32':dict(angular_max=32,points_per_interval=32,frequency_edges=EDGES,steps=24000),
    'time':dict(angular_max=16,points_per_interval=32,frequency_edges=EDGES,steps=48000),
    'range':dict(angular_max=16,points_per_interval=32,frequency_edges=EDGES[:-1],steps=24000),
    'quadrature':dict(angular_max=16,points_per_interval=48,frequency_edges=EDGES,steps=24000),
    'n32_range':dict(angular_max=32,points_per_interval=32,frequency_edges=EDGES+(512.,),steps=24000),
}


def native(x):
    if isinstance(x,np.ndarray):return native(x.tolist())
    if isinstance(x,np.generic):return native(x.item())
    if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
    if isinstance(x,dict):return {k:native(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [native(v) for v in x]
    return x


def decode(x):
    if isinstance(x,dict) and set(x)=={'real','imag'}:return complex(x['real'],x['imag'])
    if isinstance(x,dict):return {k:decode(v) for k,v in x.items()}
    if isinstance(x,list):return [decode(v) for v in x]
    return x


def tail_correct(run,coefficients):
    rows=[];tail=zeta(3,run['parameters']['angular_max']+1)
    for row,coefficient in zip(run['rows'],coefficients):
        states={}
        for name,state in row['states'].items():
            barred=dict(state['barred'])
            barred['bar_rho']+=coefficient['angular_energy_coefficient']*tail
            barred['bar_p_parallel']+=coefficient['angular_pressure_coefficient']*tail
            tensor=physical_source(barred)
            tensor['T_hatT_hatz']=state['physical']['T_hatT_hatz']
            tensor['null_plus']=tensor['radial_null']+2*tensor['T_hatT_hatz']
            tensor['null_minus']=tensor['radial_null']-2*tensor['T_hatT_hatz']
            states[name]=tensor
        rows.append({'q':row['q'],'states':states})
    return rows


def differences(first,second):
    return [{'q':a['q'],'difference':{name:{key:b['states'][name]['physical'][key]-a['states'][name]['physical'][key]
               for key in COMPONENTS} for name in a['states']}} for a,b in zip(first['rows'],second['rows'])]


def interface_checks(background):
    exact={}
    def zero(name,value):
        if isinstance(value,sp.MatrixBase):assert value.applyfunc(sp.simplify)==sp.zeros(*value.shape),name
        else:assert sp.simplify(value)==0,name
        exact[name]='0'
    f,t=sp.symbols('f t',real=True,nonnegative=True);u,v=sp.symbols('u v',real=True)
    R=u+sp.I*v;s=sp.sqrt(f*(1-f))
    CH=sp.Matrix([[f,-sp.I*s],[sp.I*s,1-f]])
    zero('affine_horizon_pure_covariance',CH**2-CH)
    M=sp.Matrix([[0,1,0],[R,0,sp.sqrt(t)]])
    C=M*sp.diag(CH,0)*M.conjugate().T
    expected=sp.Matrix([[1-f,sp.I*s*sp.conjugate(R)],[-sp.I*s*R,f*(u*u+v*v)]])
    # f is physically in [0,1], so s is real on the allowed interval.
    C=C.xreplace({sp.conjugate(s):s})
    zero('occupied_characteristic_compression',C-expected)
    zero('CAR_current_isometry',(M*M.conjugate().T-sp.eye(2)).subs(t,1-u*u-v*v))
    zero('occupied_covariance_determinant',expected.det())
    zero('occupied_covariance_trace',sp.trace(expected).subs(v*v,1-t-u*u)-(1-f*t))
    sigma1=sp.Matrix([[0,1],[1,0]])
    zero('Rindler_zero_frequency_negative_energy',expected.subs({f:sp.Rational(1,2),u:0,v:-1})-(sp.eye(2)-sigma1)/2)
    a=sp.symbols('a_EH',positive=True)
    geometric=sp.Matrix([2*a,2*a*(1-3*sp.pi),2*a*(1-3*sp.pi)])
    zero('neck_null_matches_authenticated_geometry',geometric[0]+geometric[1]-4*a*(1-3*sp.pi/2))
    rindler=[rindler_ratio_control(w,background.surface_gravity,1.) for w in (0.,.05,.2)]
    assert max(x['error'] for x in rindler)<2e-9
    controls=[]
    for angular,w in ((1,.05),(1,.2),(2,.2),(4,.5)):
        base=background.reflection(w,angular,tolerance=1e-11)
        outer=background.reflection(w,angular,outer_factor=90.,tolerance=1e-11)
        collar=background.reflection(w,angular,inner_offset=1e-12,tolerance=1e-11)
        shifted=background.reflection(w,angular,origin=.37,tolerance=1e-11)
        expected_phase=base['reflection']*np.exp(2j*w*.37)
        phase_error=abs(shifted['reflection']-expected_phase)
        original=background.interior_covariances(w,angular,base['reflection'],points=(.2,pi/2))
        alternate=background.interior_covariances(w,angular,shifted['reflection'],points=(.2,pi/2),origin=.37,delta_q=1e-12)
        covariance_error=max(np.linalg.norm(x['covariance']-y['covariance']) for x,y in zip(original,alternate))
        values={'base':{k:v for k,v in base.items() if k!='function_evaluations'},'outer_reflection_difference':outer['reflection']-base['reflection'],
                'collar_reflection_difference':collar['reflection']-base['reflection'],
                'origin_phase_error':phase_error,'independent_covariance_collar_origin_error':covariance_error,
                'propagated_covariances':original}
        assert max(abs(values['outer_reflection_difference']),abs(values['collar_reflection_difference']))<2e-8
        assert phase_error<2e-12 and covariance_error<2e-8
        for row in original:
            assert min(row['eigenvalues'])>-2e-9 and max(row['eigenvalues'])<1+2e-9
        controls.append(values)
    return {'exact':exact,'Rindler_known_solution':rindler,'actual_reflection_and_transport':controls}


def calculate(completed_development_runs=None):
    horizon,reference,_=[authenticated_record(p) for p in INPUTS]
    saved=horizon['benchmark_Unruh_per_abs_q']
    background=ParentDirac(saved['reused_horizon_rho'],saved['reused_surface_gravity'])
    runs={}
    for name,spec in SPECS.items():
        if completed_development_runs is not None and name in completed_development_runs:
            run=completed_development_runs[name]
            for key,value in spec.items():
                actual=run['parameters'][key]
                assert native(actual)==native(value),(name,key)
            assert [r['q'] for r in run['rows']]==list(POINTS)
        else:
            print('Computing new parent-source block: '+name,file=sys.stderr,flush=True)
            old=None if name in ('n16','quadrature') else runs['n32' if name=='n32_range' else 'n16']['reflection_data']
            run=parent_angular_source(background,reflection_data=old,**spec)
        runs[name]=run
    coefficients=[sixth_order_angular_coefficient(q) for q in POINTS]
    corrected={name:tail_correct(runs[name],coefficients) for name in ('n16','n32')}
    controls={name:differences(runs['n32' if name=='n32_range' else 'n16'],runs[name])
              for name in ('time','range','quadrature','n32_range')}
    thresholds={'time':5e-7,'range':5e-6,'quadrature':1e-7,'n32_range':2e-6}
    for name,rows in controls.items():
        assert max(abs(v) for r in rows for d in r['difference'].values() for v in d.values())<thresholds[name],name
    angular=[]
    for lo,hi in zip(corrected['n16'],corrected['n32']):
        d={name:{key:hi['states'][name][key]-lo['states'][name][key] for key in COMPONENTS} for name in lo['states']}
        assert max(abs(v) for state in d.values() for v in state.values())<4e-6
        angular.append({'q':hi['q'],'difference':d})
    asymptotics=[]
    for row,c in zip(runs['n32']['rows'],coefficients):
        b=row['states'][STATE]['barred']
        asymptotics.append({'q':row['q'],'sixth_order':c,
            'measured_kappa32_energy_coefficient':32**3*b['angular_energy_remainders'][-1],
            'measured_kappa32_pressure_coefficient':32**3*b['angular_pressure_remainders'][-1]})
    neck=corrected['n32'][-1]['states'][STATE]
    previous=reference['neck_source_budget']['reference_tensor']
    assert neck['null_plus']<0 and neck['null_minus']<0 and neck['rho']<0
    maximum_norm=max(row['mode_norm_defect'] for run in runs.values() for row in run['rows'])
    assert maximum_norm<2e-10
    bmin=min(np.sqrt(row['states'][STATE]['barred']['W_jets'][0]) for row in runs['n32']['rows'])
    bounds=thermal_phase_omission_bound(32,background.surface_gravity,4.,bmin)
    assert bounds['diagonal_source_bound']<1e-16
    checks=interface_checks(background)
    return native({'schema':'NSC-UNRUH-STATE-v1',
        'status':'computed characteristic parent-state candidate and canonical massless source; complete physical sourcing remains open',
        'source_hashes':hashes(SOURCES),'input_hashes':hashes(INPUTS),
        'state':{'prescription':'empty incoming parent modes plus occupied affine-horizon covariance, transported with the actual complex Dirac reflection',
                 'frame':'real current-normalized horizon spin frame; k_z=-omega; H=m sigma1+k_z/b sigma3',
                 'covariance':'[[1-f,i sqrt(f(1-f)) R*],[-i sqrt(f(1-f)) R,f |R|²]]',
                 'f':'1/(1+exp(2*pi*omega/kappa_h))','charge_pair':'C(+omega)=I-sigma3 C(-omega)* sigma3',
                 'field':'one untwisted massless four-dimensional Dirac; S2 degeneracy4*kappa',
                 'units':'hbar=c=L_throat=1; canonical renormalization mu=1; a_EH=A_EH L_throat² remains symbolic',
                 'thermal_incoming_control':'adds f |T|² to the incoming covariance; no claim of a proved global Hartle-Hawking state'},
        'method':{'stress_owner':'unchanged angular-stress adiabatic subtraction, joint finite part and full conformal anomaly',
                  'tail':'sixth-order integrated kappa^-3 remainder times Hurwitz zeta(3,N+1), checked against resolved channels and N16/N32',
                  'phase_omission':'only numerical above omega4; exact theoretical covariance retains R at every frequency',
                  'reference_recomputed':False,'raw_reflection_grids_in_record':False,
                  'published_field_policy':'all stored source rows, angular coefficients, numerical controls, source budgets and scope compared'},
        'specifications':SPECS,'runs':{k:{key:v for key,v in run.items() if key!='reflection_data'} for k,run in runs.items()},
        'tail_corrected_sources':corrected,'angular_asymptotics':asymptotics,
        'checks':{**checks,'resolution_differences':controls,'angular_tail_refinement':angular,
                  'maximum_mode_norm_defect':maximum_norm,'thermal_phase_omission_bound':bounds,
                  'uncertainty':'measured resolution changes and an asymptotic angular-tail approximation; not a rigorous bound on the full infinite tower'},
        'neck_source_budget':{'candidate_tensor':neck,'auxiliary_reference_tensor':previous,
            'candidate_minus_reference':{k:neck[k]-previous[k] for k in ('rho','p_parallel','p_sphere','radial_null')},
            'required_two_derivative_tensor':{'rho':'2*a_EH','p_parallel':'2*a_EH*(1-3*pi)',
                'p_sphere':'2*a_EH*(1-3*pi)','T_hatT_hatz':'0','null_plus_and_minus':'4*a_EH*(1-3*pi/2)'},
            'required_null_coefficient':reference['neck_source_budget']['coefficient_of_a_EH'],
            'residual_convention':'required two-derivative tensor minus computed candidate; remaining common-action terms must supply each residual',
            'density_sign_matches_positive_Einstein_requirement':False,
            'null_sign_matches':True,'pressure_anisotropy':neck['p_parallel']-neck['p_sphere'],
            'parent_Killing_power':runs['n32']['parent_power'],
            'Einstein_coefficient_fitted':False,'source_budget_closed':False},
        'remaining_owners':['global characteristic extension and Hadamard trace hypotheses on this domain',
            'same-state actual finite-Lambda causal conversion and complementary functional',
            'massive compact modes, gauge/boundary and recursive interactions from that functional',
            'independent metric, link and relative-scale equations, including flux/evolving backreaction'],
        'scope':{'inserted_scalar_source':False,'general_Unruh_or_Hadamard_theorem_claimed_new':False,
            'global_NSC_Hadamard_theorem_established':False,'full_finite_cutoff_stress_computed':False,
            'self_sourced_geometry_derived':False,'cosmological_density_Q_derived':False,
            'candidate_selected_by_required_stress_sign':False},
        'comparison':{'fields':'all published fields','float_atol':3e-9,'float_rtol':3e-8,
            'exact':'keys, strings, structures, source/input hashes, integer fields and scope','exceptions':[]}})


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite immutable evidence')
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if expected:
        compare(expected['source_hashes'],hashes(SOURCES));compare(expected['input_hashes'],hashes(INPUTS))
    result=calculate()
    if expected:
        compare(expected,result);print('Parent source: all published fields reproduced; previous scientific generators were not run.')
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
