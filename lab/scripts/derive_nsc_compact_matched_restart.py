#!/usr/bin/env python3
"""Versioned full-grid compact restart; historical seeds are controls only."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_compact_restart_scattering import generate,load_authenticated
from recursive_horizons.nsc_compact_matched_restart import assemble_versions,CompactSeedGeneration,HISTORICAL,MATCHED,VERSIONS
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain,S3
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT='results/development/nsc-compact-matched-restart.json'
SOURCES=('scripts/derive_nsc_compact_matched_restart.py','src/recursive_horizons/nsc_compact_restart_scattering.py',
         'src/recursive_horizons/nsc_compact_matched_restart.py','tests/test_nsc_compact_matched_restart.py',
         'docs/nsc-compact-matched-restart.md')
INPUTS=('results/development/nsc-mode-resolved-cauchy-state.json','results/development/charged-compact-ctp-completion.json',
        'results/development/nsc-horizon-paired-pg-map.json','src/recursive_horizons/nsc_paired_horizon_preparation.py',
        'src/recursive_horizons/nsc_compact_ctp_neck.py','src/recursive_horizons/nsc_unruh_state.py',
        'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py','src/recursive_horizons/nsc_transmitting_dirac_domain.py',
        'scripts/derive_nsc_transmitting_boundary_binding.py')


def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def conventions():
    state=json.loads((ROOT/INPUTS[0]).read_text())
    compact=json.loads((ROOT/INPUTS[1]).read_text())
    return state,compact['runs']['base']['config']


def calculate(scattering_path):
    state,config=conventions()
    scattering,metadata,receipt=load_authenticated(ROOT,scattering_path,state['channels'],config)
    arrays=assemble_versions(scattering,config,progress=lambda i,n:print(f'New versioned compact mode evolution {i}/{n}',flush=True))
    # Old covariance is opened only after both generations have been built.
    payload=state['payload']
    if sha(payload['path'])!=payload['sha256']:raise ValueError('historical seed artifact changed')
    old=np.load(ROOT/payload['path'],allow_pickle=False)
    labels=arrays['labels'];positive=labels[:,1]==1;indices=labels[positive,2].astype(int)
    if not np.array_equal(labels[positive,5],old['frequency'][indices]) or not np.array_equal(labels[positive,6],old['quadrature_weight'][indices]):
        raise ValueError('restart altered the declared frequency quadrature')
    history=arrays['covariance_'+HISTORICAL];matched=arrays['covariance_'+MATCHED]
    historical_difference=np.linalg.norm(history[positive]-old['covariance_seed'][indices],axis=(1,2))
    matched_old_difference=np.linalg.norm(matched[positive]-old['covariance_seed'][indices],axis=(1,2))
    frame_record=json.loads((ROOT/INPUTS[2]).read_text())
    shared=np.eye(2);tol=3e-11
    seam=TransmittingDiracSeamDomain(1.,1.,np.sqrt(3*np.pi/2),1.)
    trace,inverse,_=seam.trace_map(np.ones(1))
    version_checks={}
    for convention in (HISTORICAL,MATCHED):
        selected=CompactSeedGeneration.select(arrays,start_convention=convention)
        c=selected.covariance;neg=selected.negative_covariance
        eig=np.linalg.eigvalsh(np.concatenate((c,neg)))
        modes=selected.mode_map
        T=trace[0]@c@trace[0].conj().T
        restored=inverse[0]@T@inverse[0].conj().T
        signed=shared-S3@c[arrays['signed_partner_index']].conj()@S3
        version_checks[convention]={
            'version':selected.version,
            'CAR_lower_violation':max(0.,-float(eig.min())),
            'CAR_upper_violation':max(0.,float(eig.max())-1.),
            'mode_coisometry':float(np.max(np.linalg.norm(modes@modes.swapaxes(-1,-2).conj()-shared,axis=(1,2)))),
            'T_seed_feedback':float(np.max(np.linalg.norm(restored-c,axis=(1,2)))),
            'signed_source_partner':float(np.max(np.linalg.norm(neg-signed,axis=(1,2)))),
        }
    unitary_residual=float(np.max(np.linalg.norm(arrays['unitary'].swapaxes(-1,-2).conj()@arrays['unitary']-shared,axis=(1,2))))
    rows=[]
    for channel in state['channels']:
        if channel['compact_level']==0:continue
        group=labels[:,0]==channel['index'];plus=group&positive;minus=group&(~positive)
        control=labels[plus,2].astype(int)
        reference=next(x for x in frame_record['groups'] if x['channel_index']==channel['index'])
        which=np.flatnonzero(plus&(labels[:,2]==reference['stored_control_node']))[0]
        encoded=np.asarray(reference['matched_X_real_imag']);known=encoded[:,:,0]+1j*encoded[:,:,1]
        row={
            'channel_index':channel['index'],'compact_level':channel['compact_level'],'angular_level':channel['angular_level'],
            'class':'compact_only' if channel['angular_eigenvalue']==0 else 'mixed',
            'positive_X_nodes':int(plus.sum()),'positive_Y_nodes':int(minus.sum()),
            'historical_vs_locked_X':float(np.max(np.linalg.norm(history[plus]-old['covariance_seed'][control],axis=(1,2)))),
            'matched_vs_historical_X_maximum':float(np.max(np.linalg.norm(matched[plus]-old['covariance_seed'][control],axis=(1,2)))),
            'matched_vs_historical_X_minimum':float(np.min(np.linalg.norm(matched[plus]-old['covariance_seed'][control],axis=(1,2)))),
            'matched_vs_prior_independent_modal_control':float(np.linalg.norm(matched[which]-known)),
            'new_matched_control_version':VERSIONS[MATCHED],
        }
        rows.append(row)
    independent=max(r['matched_vs_prior_independent_modal_control'] for r in rows)
    maximum_checks=max([unitary_residual,float(historical_difference.max()),independent]+
        [v for row in version_checks.values() for k,v in row.items() if k!='version'])
    if maximum_checks>tol:raise ArithmeticError(json.dumps({'unitary':unitary_residual,'historical':float(historical_difference.max()),'independent':independent,'versions':version_checks},indent=2))
    declaration={
        'schema':'NSC-COMPACT-SEED-GENERATIONS-v2',
        'start_conventions':VERSIONS,
        'basis':state['basis'],
        'scope':'canonical seed-frequency mode data, NOT common-PG spatial covariance',
        'config':config,'scattering_sha256':receipt['sha256'],
        'source_hashes':{p:sha(p) for p in SOURCES},
    }
    arrays['metadata_json']=np.frombuffer(json.dumps(declaration,sort_keys=True,separators=(',',':')).encode(),dtype=np.uint8)
    binary=deterministic_npz_bytes(arrays);digest=hashlib.sha256(binary).hexdigest()
    record={
        'schema':'NSC-COMPACT-MATCHED-RESTART-v1',
        'status':'H1 PASS: explicit historical regression and new matched_delta_q compact seed control; H2/global PG state OPEN',
        'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in (*INPUTS,payload['path'])},
        'locked_inputs':frame_record['locked_inputs'],
        'scattering_payload':receipt,'scattering_provenance':metadata,
        'seed_generations_payload':{'path':f'results/development/artifacts/nsc-compact-seed-generations.{digest}.npz','sha256':digest,'bytes':len(binary)},
        'start_conventions':{
            HISTORICAL:{'version':VERSIONS[HISTORICAL],'Frobenius_argument':'sqrt(lambda^2+(m*r_h)^2)*sqrt(2*delta_q/kappa_h)/r_h','role':'historical regression only; original seed remains immutable'},
            MATCHED:{'version':VERSIONS[MATCHED],'Frobenius_argument':'sqrt(lambda^2+(m*r_h)^2)*sqrt(2*delta_q/kappa_h)','role':'new explicitly selected matched seed control; does not claim equality with historical compact X'},
        },
        'construction':{
            'source':'same affine C_H and inherited incoming occupation; complex physical scattering on required angular signs',
            'evolution':'same inherited 6000-step commutator-free SU(2) propagator shared between both initial conventions',
            'compact_groups':20,'positive_X_nodes':int(positive.sum()),'positive_Y_nodes':int((~positive).sum()),
            'signed_negative_nodes':len(labels),'angular_zero_mode_partner_not_invented':True,
            'frequencies_and_weights':'exact existing compact _frequency_grid; no new frequency or cutoff',
            'seed_covariance_used_as_input':False,'seed_covariance_used_as_control':True,
            'angular_Z':'zero from the already verified angular boundary identity; no CAR completion selected',
            'reference_and_stress_recomputed':False,
        },
        'groups':rows,
        'verification':{
            'tolerance':tol,'historical_vs_locked_X_maximum':float(historical_difference.max()),
            'matched_vs_historical_X_maximum':float(matched_old_difference.max()),
            'matched_vs_historical_X_minimum':float(matched_old_difference.min()),
            'matched_vs_prior_independent_modal_control':independent,
            'shared_propagator_unitarity':unitary_residual,'version_checks':version_checks,
            'magnetic_and_horizon_frame_certificates_reused':frame_record['verification'],
        },
        'gate':{
            'H1':'PASS','new_matched_control_surface':VERSIONS[MATCHED],
            'historical_seed_permanently_retained':True,
            'matched_values_silently_substituted_into_old_seed':False,
            'matched_checkpoint_reproduction_policy':'cold --check recomputes all modal arrays from authenticated scattering and compares all record fields and artifact bytes',
            'H2':'OPEN: full globally normalized PG modes/current resolution not supplied by a seed-frequency table',
            'C1b_M2':'OPEN','C1b_full':'OPEN','C2':'OPEN','C3':'OPEN','C4':'OPEN',
            'massive_spatial_PG_covariance':None,'transmitting_EndpointBranchJets':None,'physical_Weyl_mismatch':None,
        },
        'scope':{'parameters_refitted':False,'new_physical_terms':False,'old_seeds_modified':False,
            'LLL_or_angular_generators_rerun':False,'old_stress_reference_generators_rerun':False,
            'scattering_reused_during_modal_verification':True,'history_selected':False,
            'metric_timestep_started':False,'stress_nulls_constraints_fabricated':False,'PDF_bumped':False,'Z3':'OUT OF SCOPE'},
        'comparison':frame_record['comparison'],
    }
    return record,binary


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--prepare-scattering',action='store_true');g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    p.add_argument('--scattering',type=Path);p.add_argument('--workers',type=int,default=4);args=p.parse_args()
    if args.prepare_scattering:
        state,config=conventions()
        result=generate(ROOT,state['channels'],config,workers=args.workers,
                        progress=lambda i,n:print(f'New signed-angular scattering inputs {i}/{n}',flush=True))
        print('SCATTERING_ARTIFACT',result.relative_to(ROOT));return 0
    if args.check:
        previous=json.loads((ROOT/OUTPUT).read_text());path=ROOT/previous['scattering_payload']['path']
    else:
        if args.scattering is None:p.error('--write requires the authenticated --scattering artifact')
        path=args.scattering if args.scattering.is_absolute() else ROOT/args.scattering
    record,data=calculate(path);artifact=ROOT/record['seed_generations_payload']['path']
    if args.check:
        compare(previous,record)
        if artifact.read_bytes()!=data:raise AssertionError('versioned seed artifact differs')
        print('H1 versioned full compact restart verified; H2 remains OPEN')
    else:
        output=ROOT/OUTPUT
        if output.exists() or artifact.exists():raise FileExistsError('refusing immutable checkpoint overwrite')
        artifact.write_bytes(data);output.write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
        print(record['status']);print(record['verification'])
    return 0


if __name__=='__main__':raise SystemExit(main())
