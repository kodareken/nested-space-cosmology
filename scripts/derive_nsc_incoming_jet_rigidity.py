#!/usr/bin/env python3
"""Persist finite jet-rigidity proof; reused physical receipts are read-only."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_jet_rigidity import symbolic_proof,uniform_induction

OUTPUT='results/development/nsc-incoming-jet-rigidity.json'
SOURCES=('src/recursive_horizons/nsc_incoming_jet_rigidity.py','scripts/derive_nsc_incoming_jet_rigidity.py',
         'tests/test_nsc_incoming_jet_rigidity.py','docs/nsc-incoming-jet-rigidity.md')
PARENT='results/development/nsc-incoming-fixed-transfer.json'
INPUTS=(PARENT,'src/recursive_horizons/nsc_incoming_cauchy_jets.py')


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()


def inputs():
    parent=json.loads((ROOT/PARENT).read_text())
    if not parent['status'].startswith('PASS:') or parent['residuals']['exact_symbolic_identities']!=0:
        raise ValueError('completed physical fixed-transfer proof required')
    for section in ('source_hashes','input_hashes'):
        for path,expected in parent[section].items():
            if digest(path)!=expected:raise ValueError('reused fixed-transfer proof owner changed: '+path)
    for item in parent['reused_artifacts']:
        if digest(item['path'])!=item['sha256']:raise ValueError('reused physical P16 artifact changed')
    endpoint=parent['checked_finite_algebra']['endpoint_pair']
    if endpoint['complex_Fourier_solution']!=[{'Ahat':'0','Rhat':'0'}]:
        raise ValueError('owned complex first-normal implication is missing')
    channel=parent['actual_source_bindings']['channel']
    if channel['index']!=14 or channel['compact_mass']*channel['angular_eigenvalue']==0:
        raise ValueError('authenticated nondegenerate retained block required')
    if (parent['reused_background']['group']!=14 or parent['reused_background']['sign']!=1
            or parent['scope']['field_radial_source_old_producer_or_energy_scan_runs']!=0):
        raise ValueError('physical P16 bridge binding differs')
    return parent,endpoint['Q1']


def make_record():
    parent,q1=inputs();proof=symbolic_proof(q1);theorem=uniform_induction()
    maps=proof['endpoint_maps']['orders'];slots=proof['normal_slot_inventory']
    return {'schema':'NSC-INCOMING-JET-RIGIDITY-v1','accountable_author':'Douglas Ek',
        'status':'PASS: necessary prepared fixed-C0 tangent conditions on all twenty incoming normal slots',
        'result':{'normal_functions':['delta partial_T^k a(z)=0','delta partial_T^k r(z)=0'],
            'normal_orders':[1,2,3,4],'mixed_slots':slots['slots'],'slot_count':slots['count'],
            'endpoint_maps':maps,'physical_pair_remainder':'O(E^-6) for six inverse-gap terms',
            'sufficient_for_full_matching':False,'finite_energy_threshold':None,'uniform_large_transfer_bound':None},
        'assumptions':theorem['hypotheses'],'checked_finite_algebra':proof,'analytic_induction':theorem,
        'reused_physical_evidence':{'fixed_transfer_record':PARENT,
            'parent_physical_bridge':parent['analytic_uniform_theorem'],
            'background':parent['reused_background'],'source_bindings':parent['actual_source_bindings'],
            'physical_generators_rerun':False},
        'residuals':{'exact_symbolic_identities':0},'verification_tolerances':{'exact_symbolic_identities':0},
        'scope':{'actual_prepared_covariance_tangent':True,'mathematical_P16_frames_not_state_normalization':True,
            'N_beta_normal_jets_fixed_as_functions_through_four':True,'bare_point_germ_sufficient':False,
            'higher_normal_jets_classified':False,'complex_Fourier_amplitudes_retained':True,
            'full_matching_required_at_all_fixed_transfers':True,'finite_set_of_pairs_sufficient':False,
            'finite_amplitude_exclusion':False,'other_preparations_or_noncompact_classes_excluded':False,
            'constraint_intersection_computed':False,'full_global_constraints_or_stationarity':'OPEN',
            'broader_NSC_architecture_excluded':False,'source_law_or_boundary_prescription_changed':False,
            'new_Hadamard_gate_or_cutoff':False,'physical_field_source_radial_producer_or_energy_scan_runs':0,
            'metric_evolution_or_publication':False},
        'source_hashes':{path:digest(path) for path in SOURCES},'input_hashes':{path:digest(path) for path in INPUTS},
        'reused_artifacts':parent['reused_artifacts'],
        'reproducer':'python3 scripts/derive_nsc_incoming_jet_rigidity.py --check'}


def check():
    result=make_record()
    if json.loads((ROOT/OUTPUT).read_text())!=result:raise ValueError('jet-rigidity proof receipt differs')
    return result


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write:
        record=make_record();(ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:record=check()
    print(json.dumps({'status':record['status'],'normal_orders':record['result']['normal_orders'],
        'slot_count':record['result']['slot_count'],'physical_pair_remainder':record['result']['physical_pair_remainder'],
        'exact_scalar_residual_count':record['checked_finite_algebra']['exact_scalar_residual_count'],'residuals':record['residuals']},indent=2))


if __name__=='__main__':main()
