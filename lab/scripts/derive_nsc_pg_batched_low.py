#!/usr/bin/env python3
"""Batched real-frequency inputs for all still unresolved covariance groups."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_batched_modes import real_mode_family
from recursive_horizons.nsc_pg_batch_projection import project_mode_family
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_compact_ctp_neck import _frequency_grid
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

SOURCES=('src/recursive_horizons/nsc_pg_batched_modes.py','src/recursive_horizons/nsc_pg_batch_projection.py',
 'src/recursive_horizons/nsc_pg_massive_modes.py','src/recursive_horizons/nsc_massive_jost_modes.py',
 'src/recursive_horizons/nsc_paired_horizon_preparation.py','src/recursive_horizons/nsc_unruh_state.py',
 'src/recursive_horizons/nsc_gauge_source.py','src/recursive_horizons/nsc_lorentzian.py',
 'src/recursive_horizons/nsc_pg_lll_preparation.py','src/recursive_horizons/nsc_transmitting_dirac_domain.py',
 'src/recursive_horizons/nsc_pg_retarded_packets.py')

def one(job):
    channel,sign,config,points,upper,cache,mesh=job;g=channel['index'];suffix='-'+mesh if mesh!='coarse' else ''
    path=Path(cache)/f'group-{g}-sign-{sign}-p{points}-upper{upper:g}{suffix}.npz'
    if path.exists():return {'group':g,'sign':sign,'path':str(path),'cached':True}
    mass=channel['compact_mass'];angular=sign*channel['angular_eigenvalue']
    if mesh!='coarse':
        from numpy.polynomial.legendre import leggauss
        step=.125 if mesh=='fine' else .0625
        edges=sorted(set([float(v) for v in np.arange(0.,min(4.,upper)+step/2,step)]+[min(v,upper) for v in (mass,1.,4.,8.,12.,16.,24.,32.,upper)]))
        energies=[];weights=[]
        for left,right in zip(edges,edges[1:]):
            if right<=left or (left>=1 and right<=mass):continue
            x,q=leggauss(12 if right<=4 else (24 if mesh=='fine' else 32))
            energies.extend(left+(x+1)*(right-left)/2);weights.extend(q*(right-left)/2)
        E=np.asarray(energies);w=np.asarray(weights)
    else:
        E,w,edges=_frequency_grid(mass,upper,points);mask=~((E>1)&(E<mass));E=E[mask];w=w[mask]
    p=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],config['horizon_offset'],config['scattering_tolerance'])
    modes=real_mode_family(E,mass,angular,p)
    F=project_mode_family(E,mass,angular,modes['mode_at_zero'],modes['mode_at_three'])
    metadata={'channel':channel,'sign':sign,'points':points,'upper':upper,'edges':edges,'mesh':mesh,
              'current_residual':modes['current_residual'],'inner_unitarity_residual':modes['inner_unitarity_residual'],
              'source':'same physical horizon/infinity preparation, vectorized independent mode equations',
              'seed_covariance_used_as_input':False}
    arrays={'energy':E,'weight':w,'projection':F,'source':modes['source'],'projector':modes['projector'],
            'metadata_json':np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)}
    path.write_bytes(deterministic_npz_bytes(arrays))
    return {'group':g,'sign':sign,'path':str(path),'current':modes['current_residual'],'unitarity':modes['inner_unitarity_residual']}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--points',type=int,default=16)
    parser.add_argument('--groups',help='comma-separated group indices; defaults to all except 0,13,14');parser.add_argument('--upper',type=float,default=40.)
    parser.add_argument('--mesh',choices=('coarse','fine','finer'),default='coarse')
    args=parser.parse_args()
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    config=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    chosen=set(map(int,args.groups.split(','))) if args.groups else set(range(1,33))-{13,14}
    signature={'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},'config':config}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-batched-low-'+tag);cache.mkdir(exist_ok=True)
    jobs=[]
    for ch in state['channels']:
        if ch['index'] not in chosen:continue
        for sign in ((1,-1) if ch['angular_eigenvalue'] else (1,)):jobs.append((ch,sign,config,args.points,args.upper,str(cache),args.mesh))
    print(f'Batched physical low families {len(jobs)}; {cache}',flush=True)
    results=[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for r in pool.map(one,jobs,chunksize=1):
            results.append(r);print(json.dumps(r),flush=True)
            (cache/'progress.json').write_text(json.dumps({'signature':signature,'results':results},indent=2)+'\n')
    print(cache,flush=True)

if __name__=='__main__':main()
