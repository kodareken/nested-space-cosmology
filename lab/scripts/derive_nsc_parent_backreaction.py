#!/usr/bin/env python3
"""Bind the recorded outward power to the stationary/evolving source decision."""
import argparse
import json
from pathlib import Path
import sys

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_spinor_bridge import weyl_matrices, pauli
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT=ROOT/'results/development/parent-backreaction-gate.json'
INPUTS=('results/development/adm-neck-source-map.json',
        'results/development/warp-local-neck-source.json',
        'results/development/spherical-action.json',
        'results/development/vacuum-matched-ctp.json')
SOURCES=('scripts/derive_nsc_parent_backreaction.py','docs/nsc-parent-backreaction-gate.md',
         'src/recursive_horizons/nsc_spinor_bridge.py',
         'scripts/check_nsc_compact_boundary_action.py','scripts/check_nsc_compact_casimir.py',
         'scripts/check_nsc_vacuum_charge_matching.py')


def calculate():
    neck,local,mass,_=(authenticated_record(p) for p in INPUTS)
    power=neck['power']['stored_parent_Killing_power']
    assert power>0
    n,q=sp.symbols('N q',positive=True)
    b=sp.symbols('b_E',real=True)
    f=n*n+q*q*b*b
    shear=sp.Matrix([[1,q*q*b/f],[0,1]])
    diagonal=sp.diag(f,n*n*q*q/f)
    original=sp.Matrix([[f,q*q*b],[q*q*b,q*q]])
    assert (shear.T*diagonal*shear-original).applyfunc(sp.simplify)==sp.zeros(2)
    w=weyl_matrices();t1,_,t3=pauli()
    reflection=sp.kronecker_product(t1,w['gamma'][0])
    for eta in (-1,1):
        p=(sp.eye(8)-eta*sp.kronecker_product(t3,w['gamma5']))/2
        assert reflection*p*reflection-p==sp.zeros(8)
    beta,a=sp.symbols('beta a_EH',positive=True)
    fn,fb,fq=-8*sp.pi*a,16*sp.pi*beta*a,8*sp.pi*a*(1+2*beta*beta)
    required_power=sp.simplify(-beta*fn-(beta*beta+1)*fb+beta*fq)
    assert required_power==0
    coefficient=-power/(16*float(sp.pi))
    return {
        'schema':'NSC-PARENT-BACKREACTION-GATE-v1',
        'status':'stationary flux obstruction for the stated free analytic branch and a necessary evolving-parent mass balance',
        'source_hashes':hashes(SOURCES),'input_hashes':hashes(INPUTS),
        'reused_power':power,
        'known_local_geometric_power':local['counted_source_total']['Killing_power'],
        'exact_imports':{'Euclidean_static_shear_identity':'0',
                         'paired_reflection_preserves_compact_domain':'0 for eta=+1 and -1',
                         'static_Einstein_required_Killing_power':str(required_power)},
        'geometric_flux_result':{
            'value':'0 on the stationary symmetry-preserving analytic vacuum branch',
            'assumptions':['stationary positive Euclidean metric with noncompact radial line and no helical time identification',
                           'covariant parity-even geometric functional and reflection-compatible domain/reference',
                           'the same analytic continuation, without extra directed state or boundary data'],
            'scope':'applies to nonlocal as well as local geometric terms under these assumptions',
            'global_Lorentzian_future_horizon_reflection_assumed':False,
            'nonperturbative_horizon_continuation_existence_proved':False},
        'static_balance':{
            'required_other_sector_power':-power,
            'geometric_remainder_can_supply_it_under_stated_symmetry':False,
            'additional_independent_empty_incoming_free_channels':'nonnegative outgoing power; cannot cancel the unchanged recorded channel',
            'interacting_or_directed_recursive_state_included_in_that_bound':False,
            'full_NSC_static_no_go_claimed':False},
        'evolution_condition':{
            'reused_mass_normalization':mass['mass_and_domain'],
            'radiating_mass':'Bondi/remaining-system mass in the matched asymptotically flat regime; not total ADM energy including radiation',
            'assumptions':['positive fixed asymptotic Einstein coefficient',
                           'fixed vacuum/magnetic reference terms',
                           'spherical mean, no additional incoming power or boundary work',
                           'asymptotic retarded time u with the inherited parent normalization'],
            'physical_energy_balance':'d(M_B)/du=-P_infinity',
            'geometric_mass_balance':'dm_B/du=-G P_infinity, G=1/(16*pi*A_EH)',
            'dimensionless_initial_tangent':'a_EH * d(mhat_B)/d(uhat) = recorded coefficient',
            'aEH_times_initial_dimensionless_mass_derivative':coefficient,
            'aEH_times_rhat_times_initial_leading_g_uu_derivative':-2*coefficient,
            'rate_is_recorded_channel_contribution':True,
            'complete_total_power_known':False,
            'asymptotic_Einstein_regime_of_full_theory_established':False,
            'Einstein_coefficient_fitted':False,
            'constant_power_extrapolated_to_all_times':False},
        'decision':{'next_branch':'evolving geometry with the existing empty-incoming state',
                    'required_work':'solve the initial constraints and coupled geometry/state evolution; do not insert a mass function into the old profile and call it solved',
                    'alternative':'a directed incoming/recursive source must be derived explicitly if a different stationary state is pursued'},
        'scope':{'radiation_or_geometry_generator_rerun':False,
                 'full_diagonal_source_residuals_closed':False,
                 'finite_parent_to_child_evolution_solved':False,
                 'cosmological_density_Q_predicted':False},
        'comparison':{'fields':'all','float_atol':3e-9,'float_rtol':3e-8,
                      'exact':'structure, strings, non-float values and source/input hashes','exceptions':[]}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true');parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args();result=calculate()
    if args.check:
        compare(json.loads(args.output.read_text()),result)
        print('Backreaction gate reproduced: symmetry/domain map and all fields agree; old radiation reused.')
    else:
        with args.output.open('x') as out:
            json.dump(result,out,sort_keys=True,indent=2,allow_nan=False);out.write('\n')
        print(f'Wrote {args.output}')


if __name__=='__main__':main()
