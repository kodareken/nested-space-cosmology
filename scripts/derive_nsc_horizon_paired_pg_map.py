#!/usr/bin/env python3
"""New angular/horizon preparation and finite-collar consistency gate only."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import block_diag

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_monopole_paired_basis import verify_locked_q4_basis,paired_mode_overlap,dirac_s2_level
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap,source_covariance
from recursive_horizons.nsc_massive_signed_preparation import PAIR_SIGNED_B
from recursive_horizons.nsc_transmitting_dirac_domain import S3
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT='results/development/nsc-horizon-paired-pg-map.json'
SOURCES=('scripts/derive_nsc_horizon_paired_pg_map.py','src/recursive_horizons/nsc_monopole_paired_basis.py',
         'src/recursive_horizons/nsc_paired_horizon_preparation.py','tests/test_nsc_monopole_paired_basis.py',
         'tests/test_nsc_paired_horizon_preparation.py','docs/nsc-horizon-paired-pg-map.md')
INPUTS=('results/development/nsc-mode-resolved-cauchy-state.json','results/development/charged-ctp-neck-source.json',
        'results/development/charged-compact-ctp-completion.json','results/development/nsc-massive-signed-preparation.json',
        'results/development/nsc-pg-lll-preparation.json','src/recursive_horizons/nsc_compact_ctp_neck.py',
        'src/recursive_horizons/nsc_charged_ctp_neck.py','src/recursive_horizons/nsc_unruh_state.py',
        'src/recursive_horizons/nsc_transmitting_dirac_domain.py','src/recursive_horizons/nsc_lorentzian.py',
        'src/recursive_horizons/nsc_common_time_bulk_split.py','src/recursive_horizons/nsc_massive_signed_preparation.py',
        'docs/nsc-unruh-state.md','docs/nsc-charged-self-sourcing-route.md','scripts/derive_nsc_transmitting_boundary_binding.py')


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def pairs(x):
    x=np.asarray(x);return np.stack((x.real,x.imag),axis=-1).tolist()


def calculate():
    state,angular,compact,signed,lll=[json.loads((ROOT/p).read_text()) for p in INPUTS[:5]]
    payload=state['payload'];assert sha(payload['path'])==payload['sha256']
    arrays=np.load(ROOT/payload['path'],allow_pickle=False)
    aconfig=angular['runs']['base']['config'];cconfig=compact['runs']['base']['config']
    basis=verify_locked_q4_basis(n_max=12,n_theta=32,n_phi=32)
    angular_residuals={
        **basis['orthonormality'],**basis['ladder'],**basis['gauge'],**basis['conjugation'],
        'signed_angular_charge_map':basis['charge_map']['stated_q_positive_residual'],
        'inverse_angular_charge_map':basis['charge_map']['negative_bundle_physical_residual'],
        'locked_eta_sigma3_identification':basis['paired_lock']['max_eta1_sigma3K_lock_residual'],
    }
    angular_error=max(v for k,v in angular_residuals.items() if k not in ('n_modes','min_diagonal'))
    tolerance=3e-11
    if angular_error>tolerance:raise ArithmeticError('magnetic angular basis failed')
    owners={}
    for massless,config in ((True,aconfig),(False,cconfig)):
        owners[massless]=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],
            config['horizon_offset'],config['reflection_tolerance'] if massless else config['scattering_tolerance'],config.get('outer_floor',60.))
    rows=[]
    for channel in state['channels'][1:]:
        index=channel['index'];m=channel['compact_mass'];lam=channel['angular_eigenvalue'];n=channel['angular_level']
        level=dirac_s2_level(4,n)
        if abs(level['lam']-lam)>tolerance or level['degeneracy']!=channel['degeneracy']:
            raise ValueError('magnetic eigenspaces do not match locked spectral labels')
        h=owners[m==0]
        lo=channel['sample_offset'];hi=lo+channel['sample_count']
        node=lo+int(np.argmin(abs(arrays['frequency'][lo:hi]-.2)))
        E=float(arrays['frequency'][node])
        # Source data are computed before the stored X is inspected.
        plus=h.at_seed(E,m,lam)
        legacy=plus if m==0 else h.at_seed(E,m,lam,legacy_compact_distance=True)
        plus_negative=h.at_seed(-E,m,lam)
        matched_frame=h.frame_residual(E,m,lam)
        legacy_frame=matched_frame if m==0 else h.frame_residual(E,m,lam,legacy_compact_distance=True)
        if lam>0:
            minus=h.at_seed(E,m,-lam);minus_negative=h.at_seed(-E,m,-lam)
            overlap=paired_mode_overlap(4,n)
            modes=block_diag(plus['seed_mode_map'],minus['seed_mode_map'])
            boundary=np.kron(np.eye(2),plus['source_covariance'])
            full=modes@boundary@modes.conj().T
            quadrature=modes@np.kron(overlap,plus['source_covariance'])@modes.conj().T
            negative=block_diag(plus_negative['seed_covariance'],minus_negative['seed_covariance'])
            signed_error=float(np.linalg.norm(negative-(np.eye(4)-PAIR_SIGNED_B@full.conj()@PAIR_SIGNED_B.conj().T)))
            Y=pairs(full[2:,2:]);Z=pairs(full[:2,2:]);Ynorm=float(np.linalg.norm(full[2:,2:]-full[:2,:2]))
            angular_composition=float(np.linalg.norm(quadrature-full))
            evaluated=[plus,minus,plus_negative,minus_negative]
        else:
            full=plus['seed_covariance'];negative=plus_negative['seed_covariance']
            signed_error=float(np.linalg.norm(negative-(np.eye(2)-S3@full.conj()@S3)))
            Y=Z=None;Ynorm=None;angular_composition=0.;evaluated=[plus,plus_negative]
        observed=arrays['covariance_seed'][node]  # control only, after construction
        matched_difference=float(np.linalg.norm(plus['seed_covariance']-observed))
        legacy_difference=float(np.linalg.norm(legacy['seed_covariance']-observed))
        eigen=np.linalg.eigvalsh(full)
        residuals={
            'matched_frame_Dirac':matched_frame['differential_norm'],
            'matched_frame_offdiagonal':matched_frame['off_diagonal_maximum'],
            'source_and_T_mode_map':max(max(e['residuals'].values()) for e in evaluated),
            'angular_identity_composition':angular_composition,
            'signed_source_mode_covariance':signed_error,
            'CAR_lower_violation':max(0.,-float(eigen.min())),
            'CAR_upper_violation':max(0.,float(eigen.max())-1.),
            'legacy_seed_control':legacy_difference,
        }
        rows.append({
            'channel_index':index,'compact_mass':m,'angular_eigenvalue':lam,'angular_level':n,
            'energy':E,'stored_control_node':node,'finite_collar_delta_q':h.collar_delta_q,
            'class':'angular_only' if m==0 else ('compact_only' if lam==0 else 'mixed'),
            'residuals':residuals,
            'matched_X_real_imag':pairs(plus['seed_covariance']),
            'matched_Y_real_imag':Y,'derived_Z_real_imag':Z,
            'matched_pair_real_imag':pairs(full),'negative_pair_real_imag':pairs(negative),
            'matched_pair_eigenvalues':eigen.tolist(),'Y_minus_X_norm':Ynorm,
            'horizon_source_real_imag':pairs(plus['source_covariance']),
            'plus_seed_mode_map_real_imag':pairs(plus['seed_mode_map']),
            'plus_PG_trace_mode_map_real_imag':pairs(plus['PG_trace_mode_map']),
            'matched_X_vs_stored':matched_difference,
            'matched_X_control_status':'PASS' if matched_difference<=tolerance else 'FAIL',
            'legacy_frame_Dirac':legacy_frame['differential_norm'],
            'legacy_frame_distance_factor':legacy_frame['distance_factor'],
            'matched_vs_legacy_X_norm':float(np.linalg.norm(plus['seed_covariance']-legacy['seed_covariance'])),
            'zero_angular_partner_not_applicable':lam==0,
            'physical_spatial_PG_covariance':None,
        })
        if index%8==0:print('New paired horizon controls evaluated through group',index,flush=True)
    maxima={name:max(row['residuals'][name] for row in rows) for name in rows[0]['residuals']}
    if max(maxima.values())>tolerance:raise ArithmeticError(json.dumps(maxima,indent=2))
    angular_rows=[r for r in rows if r['class']=='angular_only'];compact_rows=[r for r in rows if r['class']!='angular_only']
    if any(r['matched_X_control_status']!='PASS' for r in angular_rows):raise ArithmeticError('angular seed reconstruction failed')
    if any(r['matched_X_control_status']!='FAIL' or r['legacy_frame_Dirac']<=1e-6 for r in compact_rows):
        raise ArithmeticError('compact normalization conflict was not isolated')
    return {
        'schema':'NSC-HORIZON-PAIRED-PG-MAP-v1',
        'status':'C1b-H OPEN: magnetic angular and matched finite modal maps evaluated; compact finite-collar checkpoint inconsistency and global spectral completion remain',
        'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in (*INPUTS,payload['path'])},
        'locked_inputs':signed['locked_inputs'],
        'magnetic_angular_basis':basis,
        'boundary_state':{
            'prescription':'existing positive affine-frequency horizon vacuum times angular spinor-bundle identity, plus inherited angular-diagonal incoming occupation',
            'matrix':'I_angular tensor (C_H(E) direct_sum n_in(E))',
            'Z_origin':'orthonormal magnetic eigenmodes and eta-conserving radial operator; not a CAR completion',
            'horizon_partner_correlations_retained':True,
            'charge_bundle_minus_q_is_not_a_flux_refit':True,
            'signed_source_relation':'C_source(-E)=I-C_source(E)*, evaluated from the charge-neutral occupation laws and the explicit conjugate-bundle basis',
        },
        'mode_map':{
            'source':'complex scattering R_eta and T_eta from the unchanged massive/angular exterior Dirac owners',
            'sewing':'M_eta=[[0,1,0],[R_eta,0,sqrt(T_eta)]]',
            'evolution':'dF/dlogdelta=-i H_logdelta F on the existing fixed Bronnikov interior',
            'PG_conversion':'same-surface T coframe/half-density map with t=tau-F(rho), Fprime=beta/A',
            'evaluated_domain':'trapped probe chart [-1,1]; seed controls at rho=0',
            'per_group_energy_policy':'one already stored frequency nearest E=0.2; no physical frequency or history selected',
            'all_retained_frequencies_evaluated':False,
            'global_PG_spectral_completeness_verified':False,
            'incoming_basis_phase':'positive sqrt(T); incoming preparation is diagonal, so its basis phase does not alter the state',
            'X_and_Y_from_seed_input':False,
        },
        'compact_collar_identification':{
            'rho_offset':'epsilon_rho=r_h^2 delta_q+O(delta_q^2)',
            'matched_Frobenius_argument':'sqrt(lambda^2+(m r_h)^2)*sqrt(2 delta_q/kappa_h)',
            'legacy_interior_argument':'matched argument divided by r_h',
            'responsible_owner':'src/recursive_horizons/nsc_compact_ctp_neck.py: compact_mode_source interior distance',
            'exterior_inverse_radius_factor_remains_correct':True,
            'both_have_same_delta_to_zero_limit':True,
            'legacy_seed_record_overwritten':False,
            'matched_candidate_adopted_as_locked_retained_state':False,
            'maximum_compact_X_change':max(r['matched_X_vs_stored'] for r in compact_rows),
            'minimum_compact_X_change':min(r['matched_X_vs_stored'] for r in compact_rows),
        },
        'verification':{'tolerance':tolerance,'maximum_angular_residual':angular_error,
            'maximum_modal_residuals':maxima,'sampled_groups':len(rows),
            'matched_X_control_pass_groups':[r['channel_index'] for r in rows if r['matched_X_control_status']=='PASS'],
            'matched_X_control_fail_groups':[r['channel_index'] for r in rows if r['matched_X_control_status']=='FAIL']},
        'groups':rows,
        'gate':{
            'magnetic_angular_normalization':'PASS','boundary_angular_Z':'PASS: fixed by declared preparation and mode orthogonality',
            'matched_modal_map':'computed at one retained frequency per group',
            'angular_only_seed_control':'PASS','compact_seed_checkpoint_consistency':'FAIL',
            'C1b_H':'OPEN','C1b_M2':'OPEN','C1b_full':'OPEN','C2':'OPEN','C3':'OPEN','C4':'OPEN',
            'remaining':'consistent versioned compact collar/state matching and globally normalized full signed PG mode resolution before spatial covariance assembly',
            'transmitting_EndpointBranchJets':None,'physical_Weyl_mismatch':None,
            'preserved_Weyl_value':signed['gate']['preserved_Weyl_value'],
            'Weyl_coefficients':signed['gate']['locked_Weyl_coefficients'],
            'Weyl_extraction_residual':signed['gate']['Weyl_extraction_residual'],
        },
        'scope':{'seed_used_only_as_control':True,'algebraic_completion_chosen_as_state':False,'LLL_map_assigned_to_massive_groups':False,
            'basis_connection_dropped':False,'parameters_refitted':False,'new_physical_terms':False,
            'old_stress_or_reference_generators_rerun':False,'new_scattering_and_modal_control_evaluations':True,
            'prior_scientific_files_modified':False,'stress_nulls_constraints_fabricated':False,
            'history_selected':False,'metric_timestep_started':False,'PDF_bumped':False,'Z3':'OUT OF SCOPE'},
        'comparison':signed['comparison'],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group()
    g.add_argument('--check',action='store_true');g.add_argument('--write',action='store_true');a=p.parse_args()
    record=calculate();path=ROOT/OUTPUT
    if a.check:
        compare(json.loads(path.read_text()),record);print('New magnetic/horizon map and compact normalization conflict verified; C1b-H OPEN')
    elif a.write:
        if path.exists():raise FileExistsError('immutable result already exists')
        path.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n');print(record['status']);print(record['verification'])
    else:print(json.dumps(record,indent=2,sort_keys=True))
    return 0


if __name__=='__main__':raise SystemExit(main())
