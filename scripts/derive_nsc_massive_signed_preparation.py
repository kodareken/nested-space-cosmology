#!/usr/bin/env python3
"""Verify massive signed maps and the mixed-channel seed-data OPEN certificate."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_massive_signed_preparation import (
    PAIR_TO_T, PAIR_SIGNED_B, signed_map_kind, map_residuals,
    exact_pairing_identities, constant_single_block_obstruction,
    mixed_pointwise_probe, paired_seed_data_witness,
)
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT='results/development/nsc-massive-signed-preparation.json'
SOURCES=('src/recursive_horizons/nsc_massive_signed_preparation.py',
         'scripts/derive_nsc_massive_signed_preparation.py',
         'tests/test_nsc_massive_signed_preparation.py','docs/nsc-massive-signed-preparation.md')
INPUTS=('results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-pg-lll-preparation.json',
        'results/development/nsc-transmitting-boundary-remainder.json',
        'results/development/nsc-general-ks-same-action-history.json',
        'results/nsc-8-chiral-boundary.json',
        'src/recursive_horizons/nsc_chiral_boundary.py','src/recursive_horizons/nsc_spinor_bridge.py',
        'src/recursive_horizons/nsc_common_time_bulk_split.py',
        'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
        'src/recursive_horizons/nsc_lorentzian.py','docs/nsc-unruh-state.md',
        'docs/nsc-chiral-boundary.md','docs/nsc-transmitting-dirac-domain.md',
        'scripts/derive_nsc_transmitting_boundary_binding.py')


def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def pairs(a):
    a=np.asarray(a);return np.stack((a.real,a.imag),axis=-1).tolist()


def calculate():
    state,lll,remainder,homogeneous=[json.loads((ROOT/f).read_text()) for f in INPUTS[:4]]
    payload=state['payload']
    if sha(payload['path'])!=payload['sha256']:raise ValueError('seed control payload changed')
    if lll['gate']['C1b_LLL']!='PASS':raise ValueError('locked LLL preparation changed')
    if not homogeneous['gate']['nonexistence_in_declared_class_proved']:raise ValueError('old homogeneous certificate changed')
    exact=exact_pairing_identities();obstruction=constant_single_block_obstruction()
    arrays=np.load(ROOT/payload['path'],allow_pickle=False)
    rows=[];mixed=[];simple=[];witness_count=0
    for channel in state['channels'][1:]:
        mass=channel['compact_mass'];angular=channel['angular_eigenvalue']
        kind=signed_map_kind(mass,angular);checks=map_residuals(mass,angular)
        row={
            'channel_index':channel['index'],'compact_mass':mass,'angular_eigenvalue':angular,
            'map_kind':kind['kind'],'radial_intertwiner_dimension':kind['dimension'],
            'current_B_real_imag':pairs(kind['current_B']),
            'seed_B_real_imag':pairs(kind['seed_B']),
            'operator_map_residuals':checks,'operator_map_status':'PASS',
            'physical_global_PG_covariance':None,
            'physical_preparation_status':'OPEN',
            'new_physical_modes_added':False,
        }
        if kind['dimension']==4:
            mixed.append(channel['index'])
            row['local_single_block_obstruction']=mixed_pointwise_probe(mass,angular)
            start=channel['sample_offset'];stop=start+channel['sample_count']
            witnesses=[paired_seed_data_witness(x) for x in arrays['covariance_seed'][start:stop]]
            witness_count+=len(witnesses)
            keys=('same_observed_marginal_residual','known_charge_paired_marginal_residual',
                  'witness_CAR_lower_violation','witness_CAR_upper_violation')
            row['information_certificate']={
                'status':'VERIFIED: recorded marginal plus CAR and signed pairing are insufficient',
                'seed_blocks_checked':len(witnesses),
                'maximum_residuals':{k:max(w[k] for w in witnesses) for k in keys},
                'minimum_unknown_negative_block_difference':min(w['unobserved_negative_same_angular_difference'] for w in witnesses),
                'positive_frequency_covariance_real_coefficients':16,
                'observed_real_coefficients':4,
                'linear_observation_kernel_dimension':12,
                'physically_free_parameter_count_claimed':False,
                'same_action_horizon_preparation_satisfied_by_witnesses':None,
                'witness_adopted_as_state':False,
            }
            # Count only after deriving the angular pair; this is bookkeeping,
            # never the equation determining its covariance or intertwiner.
            row['mode_count_bookkeeping']={
                'before_unfolding':2*channel['degeneracy']*channel['copy_count'],
                'after_unfolding':4*(channel['degeneracy']//2)*channel['copy_count'],
                'multiplicity_used_to_define_state':False,
            }
            if row['mode_count_bookkeeping']['before_unfolding']!=row['mode_count_bookkeeping']['after_unfolding']:
                raise ArithmeticError('angular unfolding changed the mode count')
        else:
            simple.append(channel['index'])
            row['information_certificate']=None
            row['missing_evaluation']='global horizon/infinity-to-PG mode reconstruction; operator symmetry alone is not a spatial state'
        rows.append(row)
    maximum_maps={k:max(row['operator_map_residuals'][k] for row in rows) for k in rows[0]['operator_map_residuals']}
    mixed_rows=[r for r in rows if r['radial_intertwiner_dimension']==4]
    maximum_witness={k:max(r['information_certificate']['maximum_residuals'][k] for r in mixed_rows)
                     for k in mixed_rows[0]['information_certificate']['maximum_residuals']}
    tol=3e-11
    if max(*maximum_maps.values(),*maximum_witness.values())>tol:raise ArithmeticError('map or witness residual failed')
    if any(r['local_single_block_obstruction']['required_basis_connection_norm']<=tol or
           r['local_single_block_obstruction']['fixed_seed_axis_potential_mismatch']<=tol for r in mixed_rows):
        raise ArithmeticError('mixed-channel obstruction not resolved')
    if simple!=list(range(1,14))+[23] or mixed!=list(range(14,23))+list(range(24,33)):
        raise ArithmeticError('retained channel classification changed')
    return {
        'schema':'NSC-MASSIVE-SIGNED-PREPARATION-v1',
        'status':'C1b-M OPEN; exact signed-spin map and mixed-channel data insufficiency certificate verified',
        'source_hashes':{p:sha(p) for p in SOURCES},
        'input_hashes':{p:sha(p) for p in (*INPUTS,payload['path'])},
        'locked_inputs':lll['locked_inputs'],
        'chosen_route':'four-component paired angular description for mixed m/lambda; simple radial involutions where one coefficient vanishes',
        'paired_operator':{
            'W':'eta1 tensor (I+i sigma2)/sqrt2',
            'W_real_imag':pairs(PAIR_TO_T),
            'H_T4':'-i I_eta tensor (v partial_rho+v_prime/2)-Nm I_eta tensor sigma1-N lambda/r eta3 tensor sigma3',
            'v':'(N/q_PG) sigma2-beta I',
            'B_T':'eta1 tensor sigma3',
            'B_T_real_imag':pairs(PAIR_SIGNED_B),
            'differential_identity':'B_T H_T4* B_T^dagger = -H_T4',
            'symbol_identity':'B_T H_T4(p)* B_T^dagger = -H_T4(-p)',
            'seed_basis_B':'eta1 tensor sigma3 after the existing T trace pullback',
            'exact_identities':exact,
            'domain':'same common PG time and transparent rho=0 seam; no reflecting spatial walls imported',
            'full_charged_4D_gauge_conjugation_claimed':False,
            'eta_means':'paired nonzero angular eigenvalue sign, not a room or compact copy',
        },
        'single_block_certificate':{
            **obstruction,
            'scope':'constant or smooth pointwise unitary signed-frequency symmetry of the same 2x2 T operator with both m and lambda/r nonzero',
            'proof':'principal term and two radii require all three Pauli anticommutators to vanish; exact kernel is zero',
            'pointwise_connection':'position-dependent frame flip carries -i v Bprime Bdagger; it cannot be dropped in the unchanged T operator',
            'full_four_component_nonexistence_claimed':False,
        },
        'state_information_certificate':{
            'observed':'X(E) in C_pair(E)=[[X,Z],[Zdagger,Y]], E>0',
            'required':'Y(E) and Z(E), or a derived symmetry/preparation law fixing them in the same paired angular/gauge basis',
            'signed_relation':'C_pair(-E)=I-B_T C_pair(E)* B_Tdagger; opposite angular signs are exchanged',
            'signed_relation_is_control_condition_not_proved_state_invariance':True,
            'witnesses':'diag(X,0) and diag(X,I), with their signed partners; algebraic controls only',
            'seed_blocks_checked':witness_count,
            'minimum_unobserved_block_difference':float(np.sqrt(2)),
            'does_not_prove':'that the action/horizon preparation cannot determine these blocks, or that NSC has no solution',
            'next_owner':'global horizon/infinity-to-paired-PG mode map and its current/gauge normalization; derive both angular marginals and preparation cross block from existing C_H and incoming occupation',
        },
        'verification':{'maximum_map_residuals':maximum_maps,'maximum_witness_residuals':maximum_witness,'tolerance':tol},
        'groups':rows,
        'gate':{
            'simple_radial_map_groups':simple,'paired_radial_map_groups':mixed,
            'massive_groups_with_PG_covariance':[],
            'C1b_M':'OPEN','C1b_full':'OPEN','C2':'OPEN','C3':'OPEN','C4':'OPEN',
            'LLL_preparation_preserved':True,
            'transmitting_EndpointBranchJets':None,'Gamma_rest_boundary_derivative':None,
            'physical_Weyl_mismatch':None,
            'locked_Weyl_coefficients':remainder['locked_local_components'],
            'Weyl_extraction_residual':remainder['locked_extraction_residual'],
            'preserved_Weyl_value':remainder['preserved_Weyl_value'],
            'homogeneous_nonexistence_preserved':True,
        },
        'scope':{'seed_used_only_as_control':True,'seed_used_to_prepare_spatial_covariance':False,
                 'LLL_map_used_for_massive_groups':False,'multiplicity_used_as_physical_map':False,
                 'reference_momentum_pairing_used_as_state':False,'new_physical_terms':False,
                 'parameters_refitted':False,'old_generators_rerun':False,'prior_scientific_files_modified':False,
                 'stress_nulls_or_constraints_fabricated':False,'history_selected':False,
                 'metric_timestep_started':False,'optimizer_started':False,'Z3':'OUT OF SCOPE',
                 'PDF_bumped':False},
        'comparison':lll['comparison'],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group()
    g.add_argument('--check',action='store_true');g.add_argument('--write',action='store_true');a=p.parse_args()
    r=calculate();path=ROOT/OUTPUT
    if a.check:
        compare(json.loads(path.read_text()),r);print('Massive signed maps and hard data-sufficiency OPEN certificate verified')
    elif a.write:
        if path.exists():raise FileExistsError('immutable result already exists')
        path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(r['status'],r['verification'])
    else:print(json.dumps(r,indent=2,sort_keys=True))
    return 0


if __name__=='__main__':raise SystemExit(main())
