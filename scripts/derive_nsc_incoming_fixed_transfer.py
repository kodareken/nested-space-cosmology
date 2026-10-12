#!/usr/bin/env python3
"""Persist physical fixed-transfer proof; exact algebra only, no propagation."""
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_fixed_transfer import symbolic_proof,uniform_proof
OUTPUT='results/development/nsc-incoming-fixed-transfer.json'
SOURCES=('src/recursive_horizons/nsc_incoming_fixed_transfer.py','scripts/derive_nsc_incoming_fixed_transfer.py',
         'tests/test_nsc_incoming_fixed_transfer.py','docs/nsc-incoming-fixed-transfer.md')
RECORDS=('results/development/nsc-incoming-endpoint-obstruction.json','results/development/nsc-incoming-vacuum-tail-bound.json')
OWNERS=tuple('src/recursive_horizons/'+name+'.py' for name in ('nsc_incoming_endpoint_obstruction','nsc_pg_high_energy',
    'nsc_incoming_vacuum_tail_bound','nsc_paired_horizon_preparation','nsc_pg_massive_modes','nsc_retarded_radial_response',
    'nsc_upstream_compensator','nsc_spatial_reference_symbol','nsc_compatible_history_geometry','nsc_ks_spacetime_variation',
    'nsc_lorentzian','nsc_transmitting_dirac_domain'))+('docs/nsc-incoming-endpoint-obstruction.md',
    'docs/nsc-incoming-fourier-matching.md','docs/nsc-retarded-compatible-preparation.md',
    'docs/nsc-incoming-vacuum-tail-bound.md','docs/nsc-pg-massive-mode-resolution.md','docs/nsc-retarded-radial-response.md')


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def inputs():
    records=[json.loads((ROOT/p).read_text()) for p in RECORDS]
    for record in records:
        for section in ('source_hashes','input_hashes'):
            for p,h in record[section].items():
                if digest(p)!=h:raise ValueError('fixed-transfer proof dependency changed: '+p)
        if 'payload' in record and digest(record['payload']['path'])!=record['payload']['sha256']:
            raise ValueError('physical vacuum-tail proof artifact changed')
    endpoint,tail=records;channel=endpoint['actual_bindings']['channel'];config=endpoint['actual_bindings']['source_config']
    if channel['index']!=14 or not channel['compact_mass']>0 or not channel['angular_eigenvalue']>0:
        raise ValueError('authenticated block14 with nonzero mass/angular label required')
    if not all(math.isfinite(config[k]) and config[k]>0 for k in ('surface_gravity','omega')):
        raise ValueError('same positive finite original source occupation parameters required')
    if endpoint['checked_finite_algebra']['endpoint']['original_limit_identity_residual']!='0':
        raise ValueError('previous physical endpoint coefficient identity failed')
    group=next(row for row in tail['groups'] if row['group']==14)
    if (group['initial_projector_difference']!=0. or group['finite_offset_archive_accuracy_certified']
        or group['initial_scope']!='declared affine-horizon limit; finite-delta input artifacts unchanged'
        or group['finite_offset_initial_projector_term_bound'] is not None):
        raise ValueError('same exact affine-horizon projector bound required, not a finite-offset archive')
    positive=next(row for row in group['per_sign'] if row['sign']==1)
    if len(positive['I16_to_I32_upper'])!=17 or not all(math.isfinite(x) and x>=0 for x in positive['I16_to_I32_upper']):
        raise ValueError('finite positive full-collar bounds required')
    return endpoint,tail,group,positive


def make_record():
    endpoint,tail,group,positive=inputs();proof=symbolic_proof();analytic=uniform_proof()
    if proof['endpoint_pair']['complex_Fourier_solution']!=[{'Ahat':'0','Rhat':'0'}]:
        raise ArithmeticError('complex pair conditions do not force both Fourier functions to vanish')
    return {'schema':'NSC-INCOMING-FIXED-TRANSFER-v1','accountable_author':'Douglas Ek',
        'status':'PASS: physical fixed-transfer first-normal necessary condition in the declared class',
        'result':{'limit01':proof['endpoint_pair']['limit01'],'limit10':proof['endpoint_pair']['limit10'],
            'coefficient_map_determinant':proof['endpoint_pair']['coefficient_map_determinant'],
            'complex_Fourier_conditions':['Ahat(omega)=0','Rhat(omega)=0'],
            'pointwise_necessary_conditions':['delta a_T(z)=0','delta r_T(z)=0'],
            'requires_matching_at_every_fixed_transfer':True,'finite_set_of_pairs_suffices':False,
            'sufficient_for_full_C0_matching':False,'physical_pair_remainder':'O(E^-4) for the constructed four-term approximation',
            'endpoint_next_order':'O(E^-3) after the displayed E^-2 coefficient',
            'finite_energy_threshold':None,'uniform_large_omega_error_bound':None},
        'assumptions':{'energy_domain':'existing unbounded signed-real-energy continuum, with E±omega/2→+infinity at fixed real omega',
            'intrinsic_fields_identically_zero_on_Sigma':['deltaN(z)','deltabeta(z)','deltaa(z)','deltar(z)'],
            'axial_class':'smooth compact first-normal tangents and retarded metric directions; Fourier injectivity uses all real omega',
            'radial_class':'fixed compact slab away from horizon; smooth energy-independent directions flat near the upstream endpoint',
            'geometry':'positive unchanged background a,r; enough bounded coefficient derivatives through q3_prime',
            'preparation':'same physical affine vacuum component and original coherent horizon/infinity source law, unchanged upstream/past',
            'nonzero_mass_angular':True},
        'checked_finite_algebra':proof,'analytic_uniform_theorem':analytic,
        'reused_background':{'record':RECORDS[1],'group':14,'sign':1,
            'background_bound_original_energy_split':group['lower'],
            'full_collar_I16_to_I32_upper':positive['I16_to_I32_upper'],
            'uniform_partial_interval_bound':'positive horizon-to-rho majorants on the slab are dominated by the existing full-collar integral',
            'needed_local_defect_derivatives':'finite background Riccati coefficient functions are smooth on compact slab, so their required rho derivatives are bounded',
            'new_coefficients_or_bound_producer_run':False},
        'actual_source_bindings':{'channel':endpoint['actual_bindings']['channel'],
            'source_config':endpoint['actual_bindings']['source_config'],'intrinsic_geometry':'rho1; a1²=3*pi/2-4, r1²=2'},
        'residuals':{'exact_symbolic_identities':0.},'verification_tolerances':{'exact_symbolic_identities':0.},
        'scope':{'actual_prepared_covariance_response':True,'formal_ad4_only_inference':False,
            'mathematical_P16_frames_not_archived_column_normalization':True,
            'coherence_or_complement_discarded':False,'source_law_or_boundary_prescription_changed':False,
            'first_normal_functions_pointwise_necessary_only':True,'higher_normal_jets_classified':False,
            'nonzero_intrinsic_directions_covered':False,'noncompact_axial_tangents_covered':False,
            'finite_amplitude_exclusion_claimed':False,'full_global_constraints_or_stationarity':'OPEN',
            'broader_NSC_architecture_excluded':False,'new_Hadamard_gate_or_cutoff':False,
            'field_radial_source_old_producer_or_energy_scan_runs':0,'metric_evolution_or_publication':False},
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in (*RECORDS,*OWNERS)},
        'reused_artifacts':[tail['payload']],
        'reproducer':'python3 scripts/derive_nsc_incoming_fixed_transfer.py --check'}


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    record=make_record();path=ROOT/OUTPUT
    if args.write:path.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=record:raise ValueError('fixed-transfer proof receipt differs')
    print(json.dumps({'status':record['status'],'result':record['result'],
        'checked_scalar_identities':record['checked_finite_algebra']['exact_scalar_residual_count'],'residuals':record['residuals']},indent=2))

if __name__=='__main__':main()
