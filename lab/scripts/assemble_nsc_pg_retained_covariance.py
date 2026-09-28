#!/usr/bin/env python3
"""Development assembly of completed low/middle/tail components only."""
import argparse
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_packet_modes import SIGNED_PACKET_MAP as S
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

def load(path):
    with np.load(path,allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}

def integrate(F,weights,source=None,projector=None):
    w=np.asarray(weights).ravel()/(2*np.pi)
    if projector is None:G=np.einsum('n,nai,nbi->ab',w,F,F.conj())
    else:G=np.einsum('n,nai,nij,nbj->ab',w,F,projector,F.conj())
    if source is None:K=np.einsum('n,nai,ij,nbj->ab',w,F,np.diag([-.5,.5,-.5]),F.conj())
    else:K=np.einsum('n,nai,nij,nbj->ab',w,F,source-.5*projector,F.conj())
    return G,K

def main():
    p=argparse.ArgumentParser();p.add_argument('--low-cache',required=True);p.add_argument('--mid-progress',default='/tmp/nsc-remaining-mid-20260915/progress.json');p.add_argument('--points',type=int,default=16);p.add_argument('--mesh',choices=('coarse','fine','finer'),default='coarse');a=p.parse_args()
    lowroot=Path(a.low_cache)
    tail_record=json.loads((ROOT/'results/development/nsc-pg-retained-tail.json').read_text());tail=load(ROOT/tail_record['payload']['path'])
    sub_record=json.loads((ROOT/'results/development/nsc-pg-retained-subgap.json').read_text());sub=load(ROOT/sub_record['payload']['path'])
    middle={}
    if Path(a.mid_progress).exists():
        for r in json.loads(Path(a.mid_progress).read_text()):
            if r['returncode']==0:
                path=Path(r['log']).read_text().strip().splitlines()[-1]
                middle[r['group'],r['sign']]=Path(path)
    positive={};missing=[]
    for i,(group,sign) in enumerate(zip(tail['group'],tail['sign'])):
        group=int(group);sign=int(sign)
        if group==14:continue # its independently prepared mixed-window inputs remain separate
        suffix='-'+a.mesh if a.mesh!='coarse' else ''
        lowpath=lowroot/f'group-{group}-sign-{sign}-p{a.points}-upper40{suffix}.npz'
        if not lowpath.exists() or (group,sign) not in middle:missing.append((group,sign));continue
        low=load(lowpath);mid=load(middle[group,sign])
        if json.loads(mid['metadata_json'].tobytes())['right']!=tail['selected_lower'][i]:raise ValueError('middle/tail window mismatch')
        G,K=integrate(low['projection'],low['weight'],low['source'],low['projector'])
        gm,km=integrate(mid['projection'],mid['weights']);G+=gm+tail['gram'][i];K+=km+tail['centered'][i]
        indices=np.flatnonzero((sub['labels'][:,0]==group)&(sub['labels'][:,1]==sign))
        if len(indices):G+=sub['gram_positive'][indices[0]]/(2*np.pi);K+=sub['centered_positive'][indices[0]]/(2*np.pi)
        positive[group,sign]=(G,K)
    rows=[];matrices={}
    for group,sign in sorted(positive):
        partner=(group,-sign) if group!=23 else (group,sign)
        if partner not in positive:continue
        G,K=positive[group,sign];Go,Ko=positive[partner];gram=G+S@Go.conj()@S;C=.5*gram+K-S@Ko.conj()@S
        eig=np.linalg.eigvalsh(C);res=float(np.linalg.norm(gram-np.eye(8),2))
        rows.append({'group':group,'sign':sign,'CAR':res,'Cmin':float(eig.min()),'Cmax':float(eig.max()),'numerical_screen':'within3e-9' if res<3e-9 else 'needs-energy-refinement'})
        matrices[f'gram_{group}_{sign}']=gram;matrices[f'covariance_{group}_{sign}']=C
    report={'scope':'development assembly; error refinement and provenance gates still apply','computed':rows,'missing_components':missing,'full_C1b':'OPEN'}
    output=Path(tempfile.gettempdir())/f'nsc-retained-covariance-screen-{a.mesh}.json';output.write_text(json.dumps(report,indent=2)+'\n')
    if matrices:(output.with_suffix('.npz')).write_bytes(deterministic_npz_bytes(matrices))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
