#!/usr/bin/env python3
"""Evaluate the matched local warp-conversion source on the actual neck.

The existing variational density and stored spectral coefficients are reused.
This computes a summand, not a truncation estimate for the full functional.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_finite_terms import _symbolic_basis
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT=ROOT/'results/development/warp-local-neck-source.json'
INPUTS=('results/development/canonical-spectral-bridge.json',
        'results/development/adm-neck-source-map.json')
SOURCES=('scripts/derive_nsc_warp_local_source.py','docs/nsc-warp-local-neck-source.md',
         'src/recursive_horizons/nsc_finite_terms.py',
         'scripts/check_nsc_compact_boundary_action.py','scripts/check_nsc_compact_casimir.py',
         'scripts/check_nsc_vacuum_charge_matching.py')


def unit_weyl_source():
    jets, _, _, _, _, _, densities=_symbolic_basis()
    extended=tuple(row+sp.symbols(name+'3:5') for row,name in zip(jets,('N','q','r')))
    derivative=lambda f:sum(sp.diff(f,row[j])*row[j+1] for row in extended for j in range(4))
    density=densities[2]
    el=[sp.diff(density,row[0])-derivative(sp.diff(density,row[1]))
        +derivative(derivative(sp.diff(density,row[2]))) for row in jets]
    x=sp.symbols('rho',real=True)
    metric_A=1+3*x+3*(1+x*x)*(sp.atan(x)-sp.pi/2)
    radius=sp.sqrt(1+x*x)
    lapse=sp.sqrt(metric_A)
    subs={z:sp.diff(f,x,j).subs(x,0) for row,f in zip(extended,(lapse,1/lapse,radius))
          for j,z in enumerate(row)}
    n,q,r=jets[0][0],jets[1][0],jets[2][0]
    static=[sp.simplify(value.subs(subs)) for value in
            (el[0]/(q*r*r),-el[1]/(n*r*r),-el[2]/(2*n*q*r))]
    # Inside A<0 the physical future time is radial, and the static
    # diagonal-frame labels exchange with these signs.
    child={'rho':-static[1],'T01':sp.S.Zero,'p_radial':-static[0],'p_sphere':static[2]}
    assert sp.simplify(child['rho']-child['p_radial']-2*child['p_sphere'])==0
    return {key:sp.simplify(value) for key,value in child.items()}


def calculate():
    bridge=authenticated_record(INPUTS[0]);neck=authenticated_record(INPUTS[1])
    coefficient=bridge['allocation']['rows'][0]['origins']['warp_covariant_conversion']['coefficients']
    v,cw,delta_a=(coefficient[key] for key in ('V_Dirac','C_Weyl_Dirac','A_Dirac'))
    unit=unit_weyl_source()
    boost=np.array(neck['coframe_map']['child_from_PG'])
    beta=neck['geometry']['beta']
    child_cw={key:float(value)*cw for key,value in unit.items()}
    child_v={'rho':v,'T01':0.,'p_radial':-v,'p_sphere':-v}
    def transform(child):
        mat=np.array([[child['rho'],child['T01']],[child['T01'],child['p_radial']]])
        out=boost.T@mat@boost
        pg={'rho':float(out[0,0]),'T01':float(out[0,1]),
            'p_radial':float(out[1,1]),'p_sphere':child['p_sphere']}
        f={'F_N':4*np.pi*pg['rho'],'F_beta':4*np.pi*pg['T01'],
           'F_q':-4*np.pi*pg['p_radial'],'F_r':-8*np.pi*pg['p_sphere']}
        flux=-beta*f['F_N']-(beta*beta+1)*f['F_beta']+beta*f['F_q']
        return {'child_tensor':child,'PG_tensor':pg,'force_densities':f,'Killing_power':float(flux)}
    parts={'vacuum':transform(child_v),'Weyl_squared':transform(child_cw)}
    total=transform({key:child_v[key]+child_cw[key] for key in child_v})
    assert abs(total['Killing_power'])<2e-13
    trace=total['PG_tensor']['rho']-total['PG_tensor']['p_radial']-2*total['PG_tensor']['p_sphere']
    assert abs(trace-4*v)<2e-14
    known={key:neck['canonical_force_densities'][key]+total['force_densities'][key]
           for key in total['force_densities']}
    known_power=-beta*known['F_N']-(beta*beta+1)*known['F_beta']+beta*known['F_q']
    assert abs(known_power-neck['power']['stored_parent_Killing_power'])<2e-13
    return {
        'schema':'NSC-WARP-LOCAL-NECK-SOURCE-v1',
        'status':'nonzero matched vacuum and Weyl-squared source contribution evaluated on the actual neck',
        'source_hashes':hashes(SOURCES),'input_hashes':hashes(INPUTS),
        'units':neck['units'],'geometry':neck['geometry'],
        'action_convention':'S_local_J=-integral sqrt(-g) [V_J+A_J R_L+C_W,J C²+...]; original metric, not the order-reduced EFT metric',
        'matched_coefficients':coefficient,
        'unit_Weyl_squared_child_source':{key:str(value) for key,value in unit.items()},
        'unit_Weyl_trace_identity':'0',
        'parts':parts,'counted_source_total':total,
        'Einstein_coefficient_increment':delta_a,
        'Einstein_bookkeeping':'A_J is included once in the complete a_EH; its -2 A_J G tensor is not also added to the source-side total',
        'canonical_plus_counted_local_J':known,
        'known_subtotal_Killing_power':float(known_power),
        'required_source_side_forces':neck['required_source_side_forces'],
        'remaining_equation':'F_unknown=F_required(a_EH)-F_canonical-F_J(vacuum+C²)',
        'scope':{'full_J_replaced_by_local_action':False,'nonlocal_remainder_bound_established':False,
                 'Euler_and_boxR_bulk_source':'0 at the interior point; actual boundary terms remain',
                 'other_sectors_R2_coefficient_set_to_zero':False,
                 'Einstein_coefficient_fitted':False,'old_source_or_spectrum_generators_rerun':False,
                 'full_source_match_or_self_sourced_geometry_claimed':False},
        'comparison':{'fields':'all','float_atol':3e-9,'float_rtol':3e-8,
                      'exact':'structure, strings, non-float values and source/input hashes','exceptions':[]}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true');parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args();result=calculate()
    if args.check:
        compare(json.loads(args.output.read_text()),result)
        print('Local warp source reproduced: all fields agree; stored coefficients and canonical tensor reused.')
    else:
        with args.output.open('x') as out:
            json.dump(result,out,sort_keys=True,indent=2,allow_nan=False);out.write('\n')
        print(f'Wrote {args.output}')


if __name__=='__main__':main()
