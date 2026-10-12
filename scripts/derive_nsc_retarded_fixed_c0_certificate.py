#!/usr/bin/env python3
"""Compose the directed projector predicate with the exact coherent sign bound.

No field transport, source quadrature or mode preparation is executed here.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

import derive_nsc_incoming_low_energy_vacuum_enclosure as V

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/development/nsc-retarded-fixed-c0-certificate.json'
SOURCES=('scripts/derive_nsc_retarded_fixed_c0_certificate.py',
         'tests/test_nsc_retarded_fixed_c0_certificate.py',
         'docs/nsc-retarded-fixed-c0-certificate.md')
INPUTS=('results/development/nsc-incoming-low-energy-vacuum-enclosure.json',
        'scripts/derive_nsc_incoming_low_energy_vacuum_enclosure.py',
        'docs/nsc-retarded-diagonal-sign-bound.md',
        'docs/nsc-incoming-fourier-matching.md',
        'results/development/nsc-retarded-compatible-response.json',
        'src/recursive_horizons/nsc_compatible_history_geometry.py',
        'src/recursive_horizons/nsc_paired_horizon_preparation.py',
        'src/recursive_horizons/nsc_pg_massive_modes.py',
        'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
        'src/recursive_horizons/nsc_common_time_bulk_split.py',
        'src/recursive_horizons/nsc_retarded_radial_response.py',
        'src/recursive_horizons/nsc_transmitting_resolvent.py',
        'src/recursive_horizons/nsc_ks_spacetime_variation.py')


def endpoint_fraction(endpoint):
    sign,mantissa,exponent,bits=map(int,endpoint)
    if sign not in (0,1) or mantissa<0 or bits<0:
        raise ValueError('finite lossless binary interval endpoint required')
    value=Q((-1 if sign else 1)*mantissa)
    return value*2**exponent if exponent>=0 else value/Q(2**(-exponent))


def binary(value):
    if not isinstance(value,(int,float)) or not math.isfinite(value):
        raise ValueError('finite frozen parameter required')
    return Q(value)


def rational(value):return {'numerator':value.numerator,'denominator':value.denominator}

def combine(interval,energy,channel,config,pulse):
    lower,upper=map(endpoint_fraction,interval)
    if lower>upper:raise ValueError('ordered directed projector enclosure required')
    if channel['index']!=14 or channel['compact_level']!=1 or channel['angular_level']!=1 or config['magnetic_flux']!=4:
        raise ValueError('owned group14 positive-angular sector required')
    E=binary(energy);kappa=binary(config['surface_gravity']);m=binary(channel['compact_mass']);ell=binary(channel['angular_eigenvalue'])
    guards={
        'selected_energy_range':Q(11,20)<E<Q(3,5),
        'positive_source_kappa_upper':0<kappa<Q(6,25),
        'positive_mass_upper':0<m<Q(8,5),
        'closed_infinity_channel':E<m,
        'positive_angular_bounds':ell>2 and ell*ell<Q(128,25),
        'same_outer_radius':pulse['normal_outer']==.03 and binary(pulse['normal_outer'])<Q(3,100),
        'same_inner_plateau':pulse['normal_inner']==.007 and binary(pulse['normal_inner'])>Q(3,500),
        'same_axial_halfwidth':pulse['axial_halfwidth']==.06 and binary(pulse['axial_halfwidth'])>Q(1,20),
        'exact_euler_lower_check':Q(27,10)**14>10**6,
    }
    if not all(guards.values()):raise ValueError('analytic sign-bound parameter guard failed')
    generator_squared=(2*Q(8,5)**2+Q(3,4)**2)/Q(4,5)**2
    occupation_exponent=2*Q(157,50)*Q(11,20)/Q(6,25)
    rotation=2*3*Q(3,100)
    commutator=2*rotation
    coherent=Q(1,500)
    margin=Q(1,2)-commutator-coherent
    insertion=Q(20,21)*Q(1,30)*Q(3,500)**2/2
    exact_upper=-margin*insertion
    if not (generator_squared<9 and occupation_exponent>14 and margin>0 and insertion>0):
        raise ArithmeticError('rational sign derivation failed')
    proved=lower>Q(1,4)
    return {'projector_lower_exact':rational(lower),'projector_upper_exact':rational(upper),
        'strict_projector_predicate':proved,'parameter_guards':guards,
        'exact_constants':{'generator_norm_squared_upper':rational(generator_squared),
            'occupation_exponent_lower':rational(occupation_exponent),
            'occupation_upper':rational(Q(1,10**6)),
            'rotation_allowance_over_I':rational(rotation),
            'vacuum_commutator_allowance_over_I':rational(commutator),
            'full_coherent_remainder_allowance_over_I':rational(coherent),
            'negative_margin_over_I':rational(margin),'positive_insertion_lower':rational(insertion)},
        'full_coherent_response_00_strict_upper':rational(exact_upper) if proved else None,
        'response_upper_outward_float':math.nextafter(float(exact_upper),math.inf) if proved else None,
        'strict_absolute_response_lower':rational(-exact_upper) if proved else None}


def signature():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (*SOURCES,*INPUTS)}


def make_record():
    vacuum=V.check()  # authenticated accounting replay; no propagation
    pilot=json.loads((ROOT/INPUTS[4]).read_text())
    if vacuum['group']!=14 or vacuum['angular_sign']!=1:
        raise ValueError('selected positive angular projector certificate required')
    result=combine(vacuum['result']['final_Re_P01']['binary_interval'],vacuum['energy'],
                   vacuum['channel'],vacuum['source_config'],pilot['control']['support'])
    if result['strict_projector_predicate']!=vacuum['result']['strict_Re_P01_above_one_quarter']:
        raise ValueError('binary predicate and projector receipt disagree')
    proved=result['strict_projector_predicate']
    return {'schema':'NSC-RETARDED-FIXED-C0-CERTIFICATE-v1','accountable_author':'Douglas Ek',
        'status':'CERTIFIED: nonzero coherent diagonal response; specified fixed-C0 tangent excluded' if proved else 'OPEN: directed projector predicate not established',
        'source_and_input_sha256':signature(),'projector_payload':vacuum['payload'],
        'group':14,'angular_sign':1,'energy':vacuum['energy'],'certificate':result,
        'physical_scope':{'object':'first derivative of full incoming covariance Fourier kernel at Eo=Ei=selected E, spin00',
            'variation':'the unchanged compact positive radius direction s*chi(s)*w(z), U=0',
            'phase_uncertainty':'all allowed subgap reflection and coherent horizon phases covered',
            'source_coherence':'retained through the exact operator-norm remainder bound',
            'source_energy_integral':'exact first-order Fourier collapse; no energy quadrature omission',
            'projector_is_a_component_not_a_replacement_state':True,
            'linearized_fixed_C0_matching':'EXCLUDED for this specified tangent' if proved else 'OPEN',
            'other_parent_extensions':'OPEN','finite_amplitude_global_matching':'OPEN',
            'global_constraint_solution':'OPEN','extended_stationarity':'OPEN','NSC_architecture':'not excluded',
            'new_field_or_radial_solves':0,'metric_timestep':False,'stress_inserted':False,
            'source_or_coupling_change':False,'additional_action_term':False,'publication':False},
        'reproducer':'python3 scripts/derive_nsc_retarded_fixed_c0_certificate.py --check'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');a=p.parse_args()
    record=make_record()
    if a.check:
        if json.loads(OUTPUT.read_text())!=record:raise ValueError('sign certificate differs from exact composition')
    else:OUTPUT.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':record['status'],'energy':record['energy'],
        'response_bound':record['certificate']['full_coherent_response_00_strict_upper'],
        'response_upper_outward_float':record['certificate']['response_upper_outward_float']},indent=2))

if __name__=='__main__':main()
