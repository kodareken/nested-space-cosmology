#!/usr/bin/env python3
"""Reuse A2 physical fields for all remaining non-subgap packet projections."""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import inspect
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_batch_projection import project_mode_family
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes


def project_family(job):
    channel,sign,source,cache=job;path=Path(cache)/f'group-{channel}-sign-{sign}.npz'
    if path.exists():return str(path)
    with np.load(source,allow_pickle=False) as a:
        labels=a['label'];m=(labels[:,0]==channel)&(labels[:,1]==sign)
        m&=(labels[:,5]<1.)|(labels[:,5]>labels[:,3])
        indices=np.flatnonzero(m);rows=labels[indices]
        if not len(indices):raise ValueError('empty declared projection family')
        mass,angular=rows[0,3:5]
        F=project_mode_family(rows[:,5],mass,angular,a['interior'][indices,1],a['exterior'][indices,1])
        data={'indices':indices,'labels':rows,'projection':F,'source_covariance':a['source_covariance'][indices],
              'source_projector':a['source_projector'][indices]}
    path.write_bytes(deterministic_npz_bytes(data));return str(path)


def main():
    record=json.loads((ROOT/'results/development/nsc-pg-massive-mode-resolution.json').read_text())
    source=ROOT/record['payload']['path']
    if hashlib.sha256(source.read_bytes()).hexdigest()!=record['payload']['sha256']:raise ValueError('A2 mode artifact changed')
    signature={'mode_sha256':record['payload']['sha256'],
       'batch_source_sha256':hashlib.sha256((ROOT/'src/recursive_horizons/nsc_pg_batch_projection.py').read_bytes()).hexdigest(),
       'worker_sha256':hashlib.sha256(inspect.getsource(project_family).encode()).hexdigest()}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-retained-projections-'+tag);cache.mkdir(exist_ok=True)
    with np.load(source,allow_pickle=False) as a:
        families=sorted(set((int(row[0]),int(row[1])) for row in a['label'] if row[0]!=13))
    # Already computed group14 scalar projections are retained as inputs.
    for sign in (1,-1):
        target=cache/f'group-14-sign-{sign}.npz'
        if target.exists():continue
        with np.load(source,allow_pickle=False) as a:
            l=a['label'];i=np.flatnonzero((l[:,0]==14)&(l[:,1]==sign)&((l[:,5]<1.)|(l[:,5]>l[:,3])))
            values=[]
            for index in i:
                old=Path('/tmp/nsc-group14-projections-4cba66cef24aa473')/f'{index:04d}.npz'
                if not old.exists():break
                with np.load(old,allow_pickle=False) as b:
                    provenance=json.loads(b['provenance_json'].tobytes())
                    if provenance['mode_payload_sha256']!=record['payload']['sha256']:raise ValueError('old group14 projection source differs')
                    values.append(b['projection'].copy())
            if len(values)==len(i):
                target.write_bytes(deterministic_npz_bytes({'indices':i,'labels':l[i],'projection':np.array(values),
                  'source_covariance':a['source_covariance'][i],'source_projector':a['source_projector'][i]}))
    print(f'Non-subgap projection families: {len(families)}; cache {cache}',flush=True)
    with ProcessPoolExecutor(max_workers=4) as pool:
        for n,path in enumerate(pool.map(project_family,((c,s,str(source),str(cache)) for c,s in families),chunksize=1),1):
            print(f'Family {n}/{len(families)}: {Path(path).name}',flush=True)
    (cache/'context.json').write_text(json.dumps(signature,indent=2)+'\n')
    print(cache,flush=True)

if __name__=='__main__':main()
