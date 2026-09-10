#!/usr/bin/env python3
"""Place the authenticated child-frame tensor into the PG ADM source equations."""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/development/adm-neck-source-map.json'
INPUTS=('results/development/unruh-state.json','results/development/adm-source-constraints.json')
SOURCES=('scripts/map_nsc_adm_neck_source.py','docs/nsc-adm-neck-source-map.md',
         'scripts/check_nsc_compact_boundary_action.py','scripts/check_nsc_compact_casimir.py',
         'scripts/check_nsc_vacuum_charge_matching.py')


def calculate():
    old=authenticated_record(INPUTS[0]);authenticated_record(INPUTS[1])
    budget=old['neck_source_budget'];t=budget['candidate_tensor']
    assert abs(t['rho_coordinate'])<1e-12 and abs(t['sphere_radius']-1)<1e-12
    b=np.sqrt(3*np.pi/2);f=b*b-1
    boost=np.array([[b,-1.],[-1.,b]])/np.sqrt(f)
    child=np.array([[t['rho'],t['T_hatT_hatz']],[t['T_hatT_hatz'],t['p_parallel']]])
    pg=boost.T@child@boost
    rho,j,p=float(pg[0,0]),float(pg[0,1]),float(pg[1,1])
    forces={'F_N':4*np.pi*rho,'F_beta':4*np.pi*j,
            'F_q':-4*np.pi*p,'F_r':-8*np.pi*t['p_sphere']}
    power=-b*forces['F_N']-(b*b+1)*forces['F_beta']+b*forces['F_q']
    frame_power=-4*np.pi*f*t['T_hatT_hatz']
    trace=rho-p-2*t['p_sphere']
    assert abs(power-budget['parent_Killing_power'])<2e-14
    assert abs(frame_power-budget['parent_Killing_power'])<2e-14
    assert abs(trace-t['trace'])<2e-14
    v,a=sp.symbols('beta a_EH',positive=True)
    transform=sp.Matrix([[v,-1],[-1,v]])/sp.sqrt(v*v-1)
    required_child=sp.diag(2*a,2*a*(1-2*v*v))
    required_pg=sp.Matrix([[-2*a,4*v*a],[4*v*a,-2*a*(1+2*v*v)]])
    residual=(transform.T*required_child*transform-required_pg).applyfunc(sp.simplify)
    assert residual==sp.zeros(2)
    return {'schema':'NSC-ADM-NECK-SOURCE-MAP-v1',
        'status':'authenticated canonical tensor transported into the four PG source equations',
        'source_hashes':hashes(SOURCES),'input_hashes':hashes(INPUTS),
        'units':old['state']['units'],
        'geometry':{'rho_coordinate':0.,'r':1.,'N':1.,'q':1.,'beta':float(b),'f_minus_A':float(f)},
        'coframe_map':{'child_time':'dT=-d rho/sqrt(f)',
                       'child_space':'sqrt(f) dz=sqrt(f) d tau+beta d rho/sqrt(f)',
                       'PG':'theta0=d tau; theta1=d rho+beta d tau',
                       'child_from_PG':boost.tolist()},
        'canonical_child_tensor':{'rho':t['rho'],'T01':t['T_hatT_hatz'],
                                  'p_radial':t['p_parallel'],'p_sphere':t['p_sphere']},
        'canonical_PG_tensor':{'rho':rho,'T01':j,'p_radial':p,'p_sphere':t['p_sphere']},
        'canonical_force_densities':forces,
        'required_Einstein_PG_tensor':{'rho':'-2*a_EH','T01':'4*beta*a_EH',
                                      'p_radial':'-2*(1+3*pi)*a_EH',
                                      'p_sphere':'2*(1-3*pi)*a_EH'},
        'required_source_side_forces':{'F_N':'-8*pi*a_EH','F_beta':'16*pi*beta*a_EH',
                                       'F_q':'8*pi*(1+3*pi)*a_EH','F_r':'16*pi*(3*pi-1)*a_EH'},
        'required_frame_identity':'0',
        'source_balance':'F_canonical+F_rest=F_required(a_EH); Einstein action contributes -F_required',
        'power':{'stored_parent_Killing_power':budget['parent_Killing_power'],
                 'from_ADM_Ward_flux':float(power),'from_child_frame':float(frame_power),
                 'Ward_flux_residual':float(power-budget['parent_Killing_power'])},
        'trace_residual':float(trace-t['trace']),
        'scope':{'Einstein_coefficient_fitted':False,'old_state_or_stress_generator_rerun':False,
                 'remaining_source_computed':False,'self_sourcing_claimed':False,
                 'parent_power_identified_with_cosmological_Q':False},
        'comparison':{'fields':'all','float_atol':3e-9,'float_rtol':3e-8,
                      'exact':'structure, strings, non-float values and source/input hashes','exceptions':[]}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true');parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args();result=calculate()
    if args.check:
        compare(json.loads(args.output.read_text()),result)
        print('ADM neck source map reproduced: all fields agree; old tensor and power reused.')
    else:
        with args.output.open('x') as out:
            json.dump(result,out,sort_keys=True,indent=2,allow_nan=False);out.write('\n')
        print(f'Wrote {args.output}')


if __name__=='__main__':main()
