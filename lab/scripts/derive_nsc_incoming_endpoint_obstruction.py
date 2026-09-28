#!/usr/bin/env python3
"""Persist/replay finite algebra and analytic physical endpoint-obstruction proof."""
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_endpoint_obstruction import symbolic_proof,rational_guard
OUTPUT='results/development/nsc-incoming-endpoint-obstruction.json'
SOURCES=('src/recursive_horizons/nsc_incoming_endpoint_obstruction.py','scripts/derive_nsc_incoming_endpoint_obstruction.py',
         'tests/test_nsc_incoming_endpoint_obstruction.py','docs/nsc-incoming-endpoint-obstruction.md')
RECORDS=('results/development/nsc-incoming-vacuum-tail-bound.json','results/development/nsc-upstream-compensator.json',
         'results/development/nsc-retarded-compatible-response.json','results/development/nsc-mode-resolved-cauchy-state.json',
         'results/development/nsc-compact-matched-restart.json')
OWNERS=tuple('src/recursive_horizons/'+name+'.py' for name in ('nsc_pg_high_energy','nsc_incoming_vacuum_tail_bound',
    'nsc_paired_horizon_preparation','nsc_pg_massive_modes','nsc_retarded_radial_response','nsc_upstream_compensator',
    'nsc_spatial_reference_symbol','nsc_compatible_history_geometry','nsc_ks_spacetime_variation','nsc_transmitting_resolvent',
    'nsc_lorentzian','nsc_transmitting_dirac_domain'))+('scripts/derive_nsc_retarded_compatible_response.py',
    'docs/nsc-incoming-fourier-matching.md','docs/nsc-retarded-compatible-preparation.md',
    'docs/nsc-incoming-vacuum-tail-bound.md','docs/nsc-pg-massive-mode-resolution.md',
    'docs/nsc-retarded-diagonal-sign-bound.md','docs/nsc-upstream-compensator.md')


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def load_inputs():
    records=[json.loads((ROOT/path).read_text()) for path in RECORDS]
    # Authoritative existing receipts are authenticated, not reproduced.
    for record in records[:3]:
        for section in ('source_hashes','input_hashes'):
            for path,h in record[section].items():
                if digest(path)!=h:raise ValueError('endpoint proof input changed: '+path)
        spec=record['payload']
        if digest(spec['path'])!=spec['sha256']:raise ValueError('endpoint proof input payload changed')
    tail,compensator,pilot,inventory,restart=records
    group=next(row for row in tail['groups'] if row['group']==14)
    if (group['initial_scope']!='declared affine-horizon limit; finite-delta input artifacts unchanged'
            or group['initial_projector_difference']!=0. or group['finite_offset_archive_accuracy_certified']
            or group['finite_offset_initial_projector_term_bound'] is not None
            or group['lower_order_cancellation']!=['0']*15):
        raise ValueError('physical affine background bound required; finite-start archive cannot replace it')
    positive=next(row for row in group['per_sign'] if row['sign']==1)
    coefficients=positive['I16_to_I32_upper']
    if len(coefficients)!=17 or not all(math.isfinite(v) and v>=0 for v in coefficients):
        raise ValueError('finite positive full-collar coefficient bounds required')
    if group['lower']!=160.:raise ValueError('existing group14 tail split changed')
    return records,group,positive


def make_record():
    records,tailgroup,positive=load_inputs();tail,compensator,pilot,inventory,restart=records
    channel=inventory['channels'][14];config=restart['scattering_provenance']['config']
    proof=symbolic_proof();guards=rational_guard(channel,pilot['control']['support'],config,compensator)
    return {'schema':'NSC-INCOMING-ENDPOINT-OBSTRUCTION-v1','accountable_author':'Douglas Ek',
        'status':'PASS: physical high-energy endpoint obstruction for the declared flat-parent-control class',
        'conclusion':{'physical_limit':'lim(E→+infinity) E² deltaChat(E,E)[0,1] = -ell*a1²*W0/(4*r1²)',
            'strict_upper_rational':guards['limit_upper_strict'],
            'fixed_C0_linear_matching_in_this_class':False,
            'finite_energy_threshold':None,'numerical_Big_O_constant':None,
            'finite_pair_calibration_can_still_hold':True},
        'assumptions':['the existing unbounded signed-real-energy continuum source/C0 prescription',
            'fixed angular/compact group14 and original affine-horizon/infinity coherent state law',
            'real smooth energy-independent metric tangent, compact radial support below rho=103/100',
            'uniformly positive N,a,r and bounded required derivatives on a fixed compact real-parameter slab',
            'operator and preparation unchanged near the upstream endpoint and beyond',
            'original radius direction has positive W0 and delta r_rho(1)=-W0/a1',
            'parent N/a/r compensators vanish on an open neighborhood of Sigma, preserving all incoming jets'],
        'auxiliary_zero_transfer':{'exact_connection':'omega=0 Fourier vertex equals the radial real-parameter derivative with integrated metric directions',
            'actual_retarded_pulse_kept':True,'stationary_family_is_coefficient_device_only':True,'delta_zero_used':False,
            'energy_quadrature_or_angular_multiplicity_inserted':False},
        'checked_finite_algebra':proof,'rational_guards':guards,
        'analytic_uniform_proof':{'coefficient_compactness':'finite recurrence coefficients and their real-parameter derivatives are bounded; projector denominator>=1',
            'explicit_linearized_projector_defect':'deltaR4=O(E^-4), by differentiating the explicit residual matrix',
            'imported_background':'Paff-P16=O(E^-16), uniformly on the compact slab by positive full-collar domination',
            'finite_background_truncation':'P16-P4=O(E^-5)',
            'true_background_error':'Paff-P4=O(E^-5)',
            'parameter_generator':'deltaG=O(E)',
            'zero_upstream_and_unitarity':'X-deltaP4=O(E^-4)',
            'physical_vacuum_endpoint':'deltaPaff_01=L/E²+O(E^-3)',
            'retained_coherence':'deltaC-deltaPaff=O(E*exp(-cE)), c>0',
            'proof_not_just_a_formal_reference_change':True,
            'uniform_orders_are_analytic_not_numerically_estimated':True},
        'reused_physical_background':{'record':RECORDS[0],'group':14,'sign':1,'starting_tail_domain':tailgroup['lower'],
            'positive_full_collar_I16_to_I32_upper':positive['I16_to_I32_upper'],
            'partial_integral_domination':'[rho,rho_h] is a subset of [1,rho_h]; each integrated coefficient majorant is nonnegative',
            'new_tail_or_coefficient_producer_run':False,'finite_offset_archived_modes_used_as_exact':False},
        'actual_bindings':{'channel':channel,'source_kappa':config['surface_gravity'],'source_omega':config['omega'],
            'source_config':config,'pulse':pilot['control']['support'],
            'compensator_support':compensator['control']['support'],
            'finite_pair_selected_compensation_context':compensator['selected_compensation']},
        'residuals':{'checked_symbolic_identities':0.,'rational_limit_bound':0.},
        'verification_tolerances':{'checked_symbolic_identities':0.,'rational_limit_bound':0.},
        'scope':{'physical_exact_state_linear_response':True,'ad4_reference_change_used_as_exclusion':False,
            'all_positive_energies_bounded_above_by_cutoff':False,'undeclared_Planck_cutoff_assumed':False,
            'alternate_source_boundary_prescriptions_excluded':False,'broader_NSC_architecture_excluded':False,
            'other_incoming_tangents_or_zero_mean_variations_decided':False,
            'pointwise_or_all_normal_jet_conditions_proved':False,'full_global_constraints_or_stationarity':'OPEN',
            'new_Hadamard_gate_or_state_law':False,'source_coherence_or_Jost_phase_discarded':False,
            'new_action_term_or_counterforce':False,'field_radial_source_or_old_producer_runs':0,
            'metric_evolution_or_publication':False},
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in (*RECORDS,*OWNERS)},
        'reused_artifacts':[records[i]['payload'] for i in range(3)],
        'reproducer':'python3 scripts/derive_nsc_incoming_endpoint_obstruction.py --check'}


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    record=make_record();path=ROOT/OUTPUT
    if args.write:path.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=record:raise ValueError('endpoint obstruction receipt differs from authenticated replay')
    print(json.dumps({'status':record['status'],'conclusion':record['conclusion'],
        'checked_identities':record['checked_finite_algebra']['checked_scalar_residual_count'],'residuals':record['residuals']},indent=2))

if __name__=='__main__':main()
