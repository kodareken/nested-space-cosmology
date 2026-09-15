#!/usr/bin/env python3
"""New real-frequency covariance inputs on unresolved smooth panels only."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import inspect
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_compact_ctp_neck import _frequency_grid
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_pg_massive_modes import MassivePGModeResolution
from recursive_horizons.nsc_pg_packet_modes import project_computed_modes
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

SOURCES=('src/recursive_horizons/nsc_paired_horizon_preparation.py','src/recursive_horizons/nsc_massive_jost_modes.py',
 'src/recursive_horizons/nsc_pg_massive_modes.py','src/recursive_horizons/nsc_pg_packet_modes.py',
 'src/recursive_horizons/nsc_compact_ctp_neck.py','src/recursive_horizons/nsc_unruh_state.py',
 'src/recursive_horizons/nsc_lorentzian.py','src/recursive_horizons/nsc_gauge_source.py',
 'src/recursive_horizons/nsc_pg_lll_preparation.py','src/recursive_horizons/nsc_common_time_bulk_split.py',
 'src/recursive_horizons/nsc_transmitting_dirac_domain.py')

def one(job):
    E,w,config,mass,angular,cache=job
    key=hashlib.sha256((float(E).hex()+','+float(w).hex()).encode()).hexdigest()[:24]
    path=Path(cache)/(key+'.npz')
    if path.exists():return str(path)
    p=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],config['horizon_offset'],config['scattering_tolerance'])
    j=solve_jost(p.background,E,mass,angular)
    mode=MassivePGModeResolution(p).section(E,mass,angular,[0.],[3.],jost=j)
    F,_=project_computed_modes(E,mass,angular,mode.interior[0],mode.exterior[0])
    values={'energy':np.array(E),'weight':np.array(w),'projection':F,'source':mode.source_covariance,'projector':mode.source_projector}
    path.write_bytes(deterministic_npz_bytes(values));return str(path)

def main():
    p=argparse.ArgumentParser();p.add_argument('--channel',type=int,required=True);p.add_argument('--sign',type=int,choices=(-1,1),default=1)
    p.add_argument('--points',type=int,default=16);p.add_argument('--upper',type=float,default=40.);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--panels',help='comma-separated existing panel bounds, e.g. 0.25:1,4:8');a=p.parse_args()
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    config=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    channel=next(c for c in state['channels'] if c['index']==a.channel)
    if a.channel in (0,13) or (channel['angular_eigenvalue']==0 and a.sign==-1):raise ValueError('choose an unresolved admissible retained family')
    mass=channel['compact_mass'];angular=a.sign*channel['angular_eigenvalue']
    E,w,edges=_frequency_grid(mass,a.upper,a.points)
    mask=~((E>1)&(E<mass));E=E[mask];w=w[mask]
    if a.panels:
        chosen=[tuple(map(float,item.split(':'))) for item in a.panels.split(',')]
        allowed=list(zip(edges,edges[1:]))
        if any(not any(abs(lo-x)<1e-12 and abs(hi-y)<1e-12 for x,y in allowed) for lo,hi in chosen):
            raise ValueError('refinement must use already declared panel endpoints')
        mask=np.zeros(len(E),bool)
        for lo,hi in chosen:mask|=(E>lo)&(E<hi)
        E=E[mask];w=w[mask]
    signature={'sources':{s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in SOURCES},
       'node_source_sha256':hashlib.sha256(inspect.getsource(one).encode()).hexdigest(),'config':config,
       'channel':channel,'angular_sign':a.sign}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-low-covariance-'+tag);cache.mkdir(exist_ok=True)
    print(f'Real panels: group {a.channel}, sign {a.sign}, {len(E)} new quadrature nodes; {cache}',flush=True)
    paths=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for n,path in enumerate(pool.map(one,((e,v,config,mass,angular,str(cache)) for e,v in zip(E,w)),chunksize=1),1):
            paths.append(path)
            if n%16==0 or n==len(E):print(f'Real nodes {n}/{len(E)}',flush=True)
    rows=[]
    for path in paths:
        with np.load(path,allow_pickle=False) as f:rows.append({k:f[k].copy() for k in f.files})
    arrays={key:np.array([r[key] for r in rows]) for key in rows[0]}
    arrays['metadata_json']=np.frombuffer(json.dumps({**signature,'points':a.points,'upper':a.upper,'edges':edges,
        'scope':'smooth real panels excluding the already-computed compact [1,m] subgap'},sort_keys=True).encode(),np.uint8)
    suffix='-'+hashlib.sha256(a.panels.encode()).hexdigest()[:8] if a.panels else ''
    output=cache/f'panels-p{a.points}-upper{a.upper:g}{suffix}.npz';output.write_bytes(deterministic_npz_bytes(arrays));print(output,flush=True)

if __name__=='__main__':main()
