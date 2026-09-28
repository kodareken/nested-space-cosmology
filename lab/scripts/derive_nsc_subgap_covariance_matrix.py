#!/usr/bin/env python3
"""Replay all inherited real subgap covariance matrices through Pauli probes."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_common_ks_trace import restrict_resolved_modes
from recursive_horizons.nsc_subgap_covariance_matrix import real_closed_source_matrix

OUTPUT=ROOT/'results/development/nsc-subgap-covariance-matrix.json'
MODE='results/development/nsc-pg-massive-mode-resolution.json'
INVENTORY='results/development/nsc-ks-source-inventory.json'
OWNED=('scripts/derive_nsc_subgap_covariance_matrix.py','docs/nsc-subgap-covariance-matrix.md',
       'src/recursive_horizons/nsc_subgap_covariance_matrix.py','src/recursive_horizons/nsc_common_ks_trace.py',
       'src/recursive_horizons/nsc_paired_horizon_preparation.py','src/recursive_horizons/nsc_subgap_history_response.py')
TOL=3e-11


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()


def payload(record):
    item=record['payload']
    if digest(item['path'])!=item['sha256']:raise ValueError('source payload changed')
    with np.load(ROOT/item['path'],allow_pickle=False) as f:return {k:f[k] for k in f.files}


def compute():
    mode_record=json.loads((ROOT/MODE).read_text());inventory=json.loads((ROOT/INVENTORY).read_text())
    mode=payload(mode_record);inputs=payload(inventory)
    config=json.loads(inputs['metadata_json'].tobytes())['config']
    labels=mode['label'];indices=np.flatnonzero((labels[:,3]>1)&(labels[:,5]>1)&(labels[:,5]<labels[:,3]))
    if len(indices)!=304:raise ValueError('all304 inherited subgap rows required')
    at_one=int(np.flatnonzero(mode['rho_interior']==1.)[0])
    energies=labels[indices,5]
    physical,_=restrict_resolved_modes(1.,energies,mode['interior'][indices,at_one])
    unsewn=mode['interior_fundamental'][indices,at_one][:,:,[1,0]]
    unsewn,_=restrict_resolved_modes(1.,energies,np.concatenate((unsewn,np.zeros((len(indices),2,1))),axis=-1))
    reflection=mode['reflection'][indices,0]
    C=mode['source_covariance'][indices]
    direct=physical@C@physical.swapaxes(-1,-2).conj()
    groups={};identity=sewing=source_error=hermiticity=0.
    for n,index in enumerate(indices):
        row=real_closed_source_matrix(energies[n],unsewn[n,:,:2],reflection[n],
            kappa=config['surface_gravity'],omega=config['omega'],mass=labels[index,3])
        result=row['covariance']
        identity=max(identity,float(np.max(abs(result-direct[n]))))
        sewn=np.column_stack((reflection[n]*unsewn[n,:,0],unsewn[n,:,1],np.zeros(2)))
        sewing=max(sewing,float(np.max(abs(sewn-physical[n]))))
        source_error=max(source_error,float(np.max(abs(row['source_covariance']-C[n]))))
        hermiticity=max(hermiticity,float(np.max(abs(result-result.conj().T))))
        group=f'{int(labels[index,0])}_{int(labels[index,1])}'
        groups[group]=groups.get(group,0)+1
    if any(v!=8 for v in groups.values()) or len(groups)!=38:raise ValueError('38 original signed subgap families required')
    passed=max(identity,sewing,source_error,hermiticity)<=TOL
    return {'schema':'NSC-SUBGAP-COVARIANCE-MATRIX-v1','accountable_author':'Douglas Ek',
        'status':('PASS' if passed else 'OPEN')+': full matrix source identity; quadrature/source accuracy OPEN',
        'rows':len(indices),'families':groups,'covariance_reconstruction_residual':identity,
        'sewing_residual':sewing,'source_covariance_residual':source_error,
        'Hermiticity_residual':hermiticity,'identity_tolerance':TOL,
        'scope':{'full_four_Pauli_components':True,'horizon_coherence_retained':True,
            'Gram_replaced_with_identity':False,'reflection_modulus_forced':False,
            'complex_bilinears_used_as_modes':False,'source_law_changed':False,
            'source_value_error_bound':None,'contour_remainder_bound':None,
            'quadrature_error_bound':None,'evolved_incoming_state_supplied':False,
            'physical_local_gate':'OPEN','new_field_or_source_runs':0,'metric_timestep':False},
        'source_hashes':{p:digest(p) for p in OWNED},
        'input_hashes':{p:digest(p) for p in (MODE,INVENTORY,mode_record['payload']['path'],inventory['payload']['path'])},
        'reproducer':'python3 scripts/derive_nsc_subgap_covariance_matrix.py --check'}


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--record',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();result=compute()
    if args.record:
        if OUTPUT.exists():raise FileExistsError('matrix identity record exists; use --check')
        OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    elif result!=json.loads(OUTPUT.read_text()):raise ValueError('matrix identity replay differs')
    print(json.dumps({k:result[k] for k in ('status','rows','covariance_reconstruction_residual','sewing_residual','source_covariance_residual','identity_tolerance')},indent=2))
