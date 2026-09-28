#!/usr/bin/env python3
"""Batched physical mode inputs on the numerical middle-energy window."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_fast_packets import FastVacuumPacketProjector
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

def main():
    p=argparse.ArgumentParser();p.add_argument('--channel',type=int,required=True);p.add_argument('--sign',type=int,choices=(-1,1),default=1)
    p.add_argument('--left',type=float,default=40.);p.add_argument('--right',type=float,default=160.);p.add_argument('--points',type=int,default=16);a=p.parse_args()
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text());ch=next(c for c in state['channels'] if c['index']==a.channel)
    config=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    if a.channel in (0,13) or (not ch['angular_eigenvalue'] and a.sign==-1):raise ValueError('choose an unresolved retained family')
    if a.left<40 and not (a.channel==14 and a.left==16):raise ValueError('lower window requires a validated boundary-mode control')
    if a.right<=a.left or (a.right-a.left)%8:raise ValueError('ordered eight-unit quadrature panels required')
    x,w=leggauss(a.points);starts=np.arange(a.left,a.right,8.)
    E=np.concatenate([v+4*(x+1) for v in starts]);weight=np.tile(4*w,len(starts))
    sources=('src/recursive_horizons/nsc_pg_fast_packets.py','src/recursive_horizons/nsc_pg_high_energy.py',
             'src/recursive_horizons/nsc_pg_lll_preparation.py','src/recursive_horizons/nsc_lorentzian.py',
             'src/recursive_horizons/nsc_transmitting_dirac_domain.py')
    metadata={'sources':{s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources},
       'config':config,'channel':ch,'angular_sign':a.sign,'left':a.left,'right':a.right,'points':a.points,
       'scope':'physical vacuum-mode approximation on this numerical interval; source corrections and boundary controls remain in acceptance ledger'}
    tag=hashlib.sha256(json.dumps(metadata,sort_keys=True).encode()).hexdigest()[:16]
    out=Path(tempfile.gettempdir())/f'nsc-mid-group{a.channel}-sign{a.sign}-{tag}.npz'
    if not out.exists():
        prep=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],config['horizon_offset'],config['scattering_tolerance'])
        F=FastVacuumPacketProjector(ch['compact_mass'],a.sign*ch['angular_eigenvalue'],prep,minimum_energy=a.left).project_many(E)
        out.write_bytes(deterministic_npz_bytes({'energies':E,'weights':weight,'projection':F,
            'metadata_json':np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)}))
    print(out,flush=True)

if __name__=='__main__':main()
