#!/usr/bin/env python3
"""Assemble the first mixed angular/compact PG covariance from both signs."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_packet_modes import SIGNED_PACKET_MAP as S
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_pg_state_covariance import equal_time_ctp_blocks

def positive_real(data,mask=None):
    E=data['energy'].ravel();mask=np.ones(len(E),bool) if mask is None else mask
    F=data['projection'][mask];w=data['weight'].ravel()[mask]/(2*np.pi)
    G=np.einsum('n,nai,nij,nbj->ab',w,F,data['projector'][mask],F.conj())
    K=np.einsum('n,nai,nij,nbj->ab',w,F,data['source'][mask]-.5*data['projector'][mask],F.conj())
    return G,K

def positive_mid(data):
    F=data['projection'];w=data['weights']/(2*np.pi);D=np.diag([-.5,.5,-.5])
    return np.einsum('n,nai,nbi->ab',w,F,F.conj()),np.einsum('n,nai,ij,nbj->ab',w,F,D,F.conj())

def load(path):
    with np.load(path,allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}

def assemble(inputs):
    positives={};low_changes={};mid_changes={};components={}
    for sign in (1,-1):
        base=inputs[f'low16_{sign}'];ref=inputs[f'low32_{sign}'];E=base['energy'].ravel()
        mask=((E>.25)&(E<1))|((E>np.pi/2)&(E<8))
        G,K=positive_real(base);oldG,oldK=positive_real(base,mask);newG,newK=positive_real(ref)
        G+=newG-oldG;K+=newK-oldK;low_changes[sign]=(newG-oldG,newK-oldK)
        mG,mK=positive_mid(inputs[f'mid24_{sign}']);cG,cK=positive_mid(inputs[f'mid16_{sign}'])
        mid_changes[sign]=(mG-cG,mK-cK)
        G+=mG+inputs[f'subgap_gram_{sign}']/(2*np.pi)+inputs[f'tail_gram_{sign}']
        K+=mK+inputs[f'subgap_centered_{sign}']/(2*np.pi)+inputs[f'tail_centered_{sign}']
        positives[sign]=(G,K)
    results={};checks={}
    for sign in (1,-1):
        G,K=positives[sign];Go,Ko=positives[-sign]
        gram=G+S@Go.conj()@S;centered=K-S@Ko.conj()@S;C=.5*gram+centered
        results[sign]={'covariance':C,'car_gram':gram,**equal_time_ctp_blocks(C,gram)}
        checks[sign]={'CAR':float(np.linalg.norm(gram-np.eye(8),2)),
          'Hermiticity':float(np.linalg.norm(C-C.conj().T,2)),
          'eigenvalue_min':float(np.linalg.eigvalsh(C).min()),'eigenvalue_max':float(np.linalg.eigvalsh(C).max()),
          'probe_bulk_correlation':float(np.linalg.norm(C[:4,4:],2)),
          'low_refinement_C':float(np.linalg.norm(.5*(low_changes[sign][0]+S@low_changes[-sign][0].conj()@S)+low_changes[sign][1]-S@low_changes[-sign][1].conj()@S,2)),
          'mid_refinement_C':float(np.linalg.norm(.5*(mid_changes[sign][0]+S@mid_changes[-sign][0].conj()@S)+mid_changes[sign][1]-S@mid_changes[-sign][1].conj()@S,2))}
    checks['signed_covariance_relation']=float(np.linalg.norm(results[-1]['covariance']-S@(results[1]['car_gram']-results[1]['covariance']).conj()@S,2))
    return results,checks

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);a=p.parse_args()
    output=Path(a.input);arrays=load(output)
    if hashlib.sha256(output.read_bytes()).hexdigest() not in output.name:raise ValueError('input digest differs')
    inputs={}
    for key,value in arrays.items():
        if '/' in key:
            group,name=key.split('/',1);inputs.setdefault(group,{})[name]=value
        else:inputs[key]=value
    results,checks=assemble(inputs)
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
