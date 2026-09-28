#!/usr/bin/env python3
"""Check the finite-endpoint/state matching contract; reuse old source records."""
import argparse
import json
from math import exp,log10,pi
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import erfc
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_state_regulator import (
    canonical_thermal,endpoint_images,matsubara_source,canonical_integral,mode_energy_derivative,
)
from recursive_horizons.nsc_influence import thermal_covariance,ground_covariance,mean
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare,native
from check_nsc_vacuum_charge_matching import hashes

OUTPUT=ROOT/'results/development/state-regulator.json'
SOURCES=('scripts/check_nsc_state_regulator.py','src/recursive_horizons/nsc_state_regulator.py',
         'docs/nsc-state-regulator.md','scripts/check_nsc_compact_boundary_action.py',
         'scripts/check_nsc_compact_casimir.py','scripts/check_nsc_vacuum_charge_matching.py')
INPUTS=('results/development/unruh-state.json','results/development/compact-matching.json',
        'results/nsc-10-influence.json','results/development/horizon-source.json')
FIELDS=('free_energy','rho','pressure','entropy','trace')


def calculate():
    parent,compact,_,horizon=[authenticated_record(p) for p in INPUTS]
    exact={}
    def zero(name,value):
        assert sp.simplify(value)==0,name
        exact[name]='0'
    T,nu,n,z=sp.symbols('T nu n z',positive=True)
    x=n*n*nu*nu/(4*T*T)
    endpoint=4*T**4/(sp.pi**2*n**4)*(1+x)*sp.exp(-x)
    density=endpoint-T*sp.diff(endpoint,T)
    zero('one_image_lapse_variation',density+3*endpoint+8*T**4*x*x*sp.exp(-x)/(sp.pi**2*n**4))
    zero('one_image_scale_trace',density+3*endpoint-nu*sp.diff(endpoint,nu))
    zero('canonical_thermal_trace',(-7*sp.pi**2*T**4/180)-T*sp.diff(-7*sp.pi**2*T**4/180,T)-7*sp.pi**2*T**4/60)
    zero('unit_residue_of_finite_proper_time_inverse',sp.limit(z*sp.exp(-z/nu**2)/z,z,0)-1)
    beta,V,c=sp.symbols('beta V c',positive=True)
    local_density=(c*beta*V)/(beta*V)
    zero('local_vacuum_shift_cancels_between_states',local_density-local_density.subs(beta,2*beta))
    # The UV divergent heat coefficients cancel in state differences on the
    # same background; the thermal winding remainder is finite, not zero.
    rows=[]
    for temperature in (.1,.25,.5,1.):
        image=endpoint_images(temperature);direct=matsubara_source(temperature)
        errors={k:direct[k]-image['finite_endpoint'][k] for k in FIELDS}
        assert max(abs(v) for v in errors.values())<2e-11
        integral=canonical_integral(temperature)
        assert abs(integral['rho']-image['canonical']['rho'])<2e-12
        step=temperature*1e-4
        F=lambda t:endpoint_images(t)['finite_endpoint']['free_energy']
        derivative=(F(temperature-2*step)-8*F(temperature-step)+8*F(temperature+step)-F(temperature+2*step))/(12*step)
        observed=F(temperature)-temperature*derivative
        variation_error=observed-image['finite_endpoint']['rho']
        assert abs(variation_error)<2e-10
        rows.append({'images':image,'direct_Matsubara':direct,'representation_difference':errors,
                     'canonical_occupation_integral':integral,'independent_density_derivative':observed,
                     'density_derivative_error':variation_error})
    assert rows[2]['images']['finite_endpoint']['rho']<0<rows[2]['images']['canonical']['rho']
    crossing=brentq(lambda t:endpoint_images(t)['finite_endpoint']['rho'],.25,1.,xtol=1e-13)

    modes=[]
    for temperature in (.25,.5,1.):
        energy=.7;row=mode_energy_derivative(energy,temperature,images=128)
        full_derivative=-4*energy*temperature*sum(exp(-((2*j+1)**2*pi*pi*temperature**2+energy**2))
                            /((2*j+1)**2*pi*pi*temperature**2+energy**2) for j in range(64))
        direct_thermal=full_derivative+erfc(energy)
        assert abs(direct_thermal-row['finite_thermal_derivative'])<2e-12
        h=np.diag([energy,-energy]);covariance=thermal_covariance(h,1/temperature)
        _,_,ground=ground_covariance(h)
        excitation=float(mean(covariance-ground,h).real)
        expected=energy*row['canonical_thermal_derivative']
        assert abs(excitation-expected)<2e-13 and excitation>0
        modes.append({'energy_over_nu':energy,'temperature_over_nu':temperature,**row,
                      'direct_Matsubara_thermal_derivative':direct_thermal,
                      'derivative_difference':direct_thermal-row['finite_thermal_derivative'],
                      'canonical_CAR_excitation_energy':excitation,'expected_excitation_energy':expected})

    temperature=horizon['benchmark_Unruh_per_abs_q']['Hawking_temperature']
    cold=endpoint_images(temperature)
    log_difference=log10(abs(cold['endpoint_minus_canonical']['rho']))
    assert log_difference<-70
    moments=compact['refined_moments']['Q2']['weight']
    full_vacuum=moments/(8*pi*pi)
    rebuilt=[]
    for row in compact['matched_coefficients']:
        value=row['V_Dirac']+row['matching_cutoff']**4/(16*pi*pi)
        assert abs(value-full_vacuum)<3e-12
        rebuilt.append({'nu':row['matching_cutoff'],'full_raw_vacuum_density':value,'difference_from_refined_moment':value-full_vacuum})
    return native({'schema':'NSC-STATE-REGULATOR-v1',
        'status':'computed state-dependent endpoint conversion and a canonical thermal-source incompatibility of the raw finite determinant',
        'source_hashes':hashes(SOURCES),'input_hashes':hashes(INPUTS),'exact_checks':exact,
        'domain':{'metric':'flat R3 times antiperiodic Euclidean thermal S1; no identification with the expanding child',
                  'field':'one massless four-component Dirac field','units':'hbar=c=k_B=1; tabulated nu=1',
                  'nu':'auxiliary light-sector matching cutoff, distinct from carrier Lambda and temperature',
                  'physical_comparison':'canonical Gaussian/CAR state from the existing influence functional'},
        'formulae':{'thermal_endpoint_free_energy':'4*T^4/pi² sum_n>=1 (-1)^n/n^4 [1-(1+x_n)exp(-x_n)]',
                   'x_n':'n²*nu²/(4*T²)',
                   'canonical_free_energy':'-7*pi²*T^4/180',
                   'thermal_C_nu_mu':'-4*T^4/pi² sum_n>=1 (-1)^n/n^4 (1+x_n)exp(-x_n)',
                   'physical_density':'F-T*dF/dT','thermal_trace':'nu^4/(2*pi²) sum_n>=1 (-1)^n exp(-x_n)',
                   'canonical_pair_excitation':'2*E/(exp(E/T)+1)>=0'},
        'thermal_source_rows':rows,'mode_residue_and_state_controls':modes,
        'negative_thermal_energy_control':{'thermal_energy_zero_T_over_nu':crossing,'root_residual':endpoint_images(crossing)['finite_endpoint']['rho'],
            'physical_temperature_limit_claimed':False,'high_T_light_limit_rho_over_nu4':-1/(16*pi*pi)},
        'cold_flat_control':{'temperature_reused':temperature,'calculation':cold,
             'log10_abs_density_conversion':log_difference,'actual_curved_source_error_bound':False},
        'raw_carrier_high_temperature_constraint':{'reused_Q2_h':moments,'reused_Lambda':2.,'reused_interval':2.,
             'full_vacuum_density':full_vacuum,'matching_split_reconstruction':rebuilt,
             'conditional_high_T_thermal_density_limit':-full_vacuum,
             'conditions':'positive finite free Euclidean modulus; no fermionic Matsubara zero mode; F_raw and T*dF_raw/dT tend to zero',
             'weight_control':'the existing positive weighted compact heat operator has exponentially suppressed large-frequency modes',
             'finite_temperature_warped_integral_rerun':False,'validity_of_physical_NSC_at_arbitrary_T_assumed':False},
        'source_matching_contract':{'allocation':'Gamma5_CTP=Gamma_light_ren_CTP[C]+Gamma_H_CTP[C]+C_nu_mu_CTP[C]',
             'C_nu_mu_can_be_declared_state_independent':False,
             'unit_pole_residue_proves_thermal_matching':False,
             'state_independent_local_terms_repair_fixed_background_thermal_difference':False,
             'canonical_completion_condition':'Delta Gamma_completion=Delta Gamma_canonical-Delta Gamma_raw for a retained canonical sector and common state/contour',
             'completion_is_a_fitted_stress':False,
             'nonthermal_curved_completion_evaluated':False,
             'normalized_relative_functional_alone_proves_full_positivity':False},
        'parent_source_reused':{'canonical_neck_tensor':parent['neck_source_budget']['candidate_tensor'],
             'parent_Killing_power':parent['neck_source_budget']['parent_Killing_power'],
             'canonical_source_unchanged':True,'source_or_scattering_generator_rerun':False,
             'full_finite_cutoff_source_identified_with_this_tensor':False},
        'remaining_owner':'common closed-time-path and initial-state prescription, with the derived state conversion and remaining sectors, before the full metric source solve',
        'scope':{'new_general_thermal_theorem_claimed':False,'new_scalar_source_inserted':False,
             'old_physical_record_modified':False,'complete_theory_disproved':False,'self_sourcing_derived':False,
             'raw_finite_modulus_alone_is_the_unchanged_canonical_state_functional':False},
        'comparison':{'fields':'all','float_atol':3e-9,'float_rtol':3e-8,
             'exact':'keys, strings, structures, hashes and scope','tiny_field_control':'log10_abs_density_conversion preserves the tiny cold mismatch','exceptions':[]}})


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite evidence')
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if expected:
        compare(expected['source_hashes'],hashes(SOURCES));compare(expected['input_hashes'],hashes(INPUTS))
    result=calculate()
    if expected:
        compare(expected,result);print('State/regulator conversion: all fields reproduced; previous scientific generators were not run.')
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
