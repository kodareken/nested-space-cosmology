#!/usr/bin/env python3
"""Physical closed-channel PG packet integral on a retarded frequency contour.

Only the unresolved subgap panel is evaluated. Other spectral intervals and
the high-energy tail are not replaced or silently assigned a covariance.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import inspect
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_complex_horizon_modes import complex_reflection,finite_pg_horizon_basis
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_pg_packet_modes import project_horizon_bases,SIGNED_PACKET_MAP
from recursive_horizons.nsc_pg_threshold_projection import contour_integral
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

DEPENDENCIES=('src/recursive_horizons/nsc_complex_horizon_modes.py',
              'src/recursive_horizons/nsc_paired_horizon_preparation.py',
              'src/recursive_horizons/nsc_pg_packet_modes.py',
              'src/recursive_horizons/nsc_pg_threshold_projection.py',
              'src/recursive_horizons/nsc_pg_massive_modes.py',
              'src/recursive_horizons/nsc_common_time_bulk_split.py',
              'src/recursive_horizons/nsc_lorentzian.py',
              'src/recursive_horizons/nsc_unruh_state.py',
              'src/recursive_horizons/nsc_gauge_source.py',
              'src/recursive_horizons/nsc_pg_lll_preparation.py',
              'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
              'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py')

def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def node_key(z):return hashlib.sha256((float(z.real).hex()+','+float(z.imag).hex()).encode()).hexdigest()[:24]

def evaluate_node(job):
    z,config,mass,angular,cache=job;path=Path(cache)/(node_key(z)+'.npz')
    if path.exists():return
    p=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],config['horizon_offset'],config['scattering_tolerance'])
    basis=finite_pg_horizon_basis(z,mass,angular,p)
    projected,_=project_horizon_bases(z,mass,angular,basis.interior,basis.exterior,rtol=2e-13,atol=2e-15)
    if z.imag:
        dual,_=project_horizon_bases(z.conjugate(),mass,angular,basis.interior_conjugate,basis.exterior_conjugate,rtol=2e-13,atol=2e-15)
    else:dual=projected.copy()
    arrays={'frequency':np.array(z),'projection':projected,'projection_conjugate':dual}
    residuals={'basis':basis.residuals}
    if z.imag>0:
        reflection=complex_reflection(z,mass,angular,p.background,radial_collar=config['horizon_offset'])
        arrays['reflection']=np.array(reflection.reflection)
        residuals['reflection']=reflection.residuals
    arrays['residuals_json']=np.frombuffer(json.dumps(residuals,sort_keys=True).encode(),np.uint8)
    path.write_bytes(deterministic_npz_bytes(arrays))

def nodes_for(left,right,kappa,points,height_fraction):
    height=kappa/height_fraction;x,_=leggauss(points);points_z=[]
    for a,b in ((left,right),(complex(left),left+1j*height),
                (left+1j*height,right+1j*height),(right+1j*height,complex(right))):
        points_z.extend(complex(a+(t+1)*(b-a)/2) for t in x)
    return points_z

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--channel',type=int,default=13)
    parser.add_argument('--sign',type=int,choices=(-1,1),default=1)
    parser.add_argument('--points',type=int,default=16)
    parser.add_argument('--height-fraction',type=float,default=4.)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    previous=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())
    config=previous['scattering_provenance']['config']
    channel=next(c for c in state['channels'] if c['index']==args.channel)
    mass=float(channel['compact_mass']);angular=args.sign*float(channel['angular_eigenvalue'])
    if mass<=1 or (angular==0 and args.sign!=1):raise ValueError('existing compact subgap panel and admissible angular sign required')
    if args.height_fraction<=2:raise ValueError('contour must remain below first source/Frobenius pole')
    signature={'source_hashes':{p:sha(p) for p in DEPENDENCIES},
               'node_function_sha256':hashlib.sha256(inspect.getsource(evaluate_node).encode()).hexdigest(),
               'config':config,'channel':channel,'angular_sign':args.sign}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-threshold-contour-'+tag);cache.mkdir(exist_ok=True)
    nodes=nodes_for(1.,mass,config['surface_gravity'],args.points,args.height_fraction)
    unique=list(dict.fromkeys(nodes))
    print(f'New subgap contour: channel {args.channel}, sign {args.sign}, {len(unique)} nodes; cache {cache}',flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,_ in enumerate(pool.map(evaluate_node,((z,config,mass,angular,str(cache)) for z in unique),chunksize=1),1):
            if i%16==0 or i==len(unique):print(f'Contour fields {i}/{len(unique)}',flush=True)
    memo={};reflections={}
    for z in unique:
        with np.load(cache/(node_key(z)+'.npz'),allow_pickle=False) as a:
            if a['frequency'].item()!=z:raise ValueError('cached frequency differs')
            memo[z]=a['projection'].copy();memo[z.conjugate()]=a['projection_conjugate'].copy()
            if z.imag>0:reflections[z]=a['reflection'].item()
    result=contour_integral(1.,mass,config['surface_gravity'],lambda z:memo[z],lambda z:reflections[z],
                            points=args.points,height=config['surface_gravity']/args.height_fraction)
    data={'gram_positive':result.gram_positive_energy,'centered_positive':result.centered_positive_energy,
          'analytic':result.analytic_integral,'smooth':result.smooth_integral}
    if angular==0:
        data['gram_signed']=(result.gram_positive_energy+SIGNED_PACKET_MAP@result.gram_positive_energy.conj()@SIGNED_PACKET_MAP)/(2*np.pi)
        data['centered_signed']=(result.centered_positive_energy-SIGNED_PACKET_MAP@result.centered_positive_energy.conj()@SIGNED_PACKET_MAP)/(2*np.pi)
    metadata={**signature,'points':args.points,'height_fraction':args.height_fraction,
              'interval':[1.,mass],'contour_nodes_real_imag':[[z.real,z.imag] for z in unique],
              'scope':'subgap panel only; no full covariance or spectral tail assigned'}
    data['metadata_json']=np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)
    output=cache/f'panel-p{args.points}-h{args.height_fraction:g}.npz';output.write_bytes(deterministic_npz_bytes(data))
    print(json.dumps({'panel':str(output),'Gram_positive_trace':float(np.trace(result.gram_positive_energy).real),
                      'Gram_positive_eigenvalues':np.linalg.eigvalsh(result.gram_positive_energy).tolist(),
                      'Hermitian_residual':float(np.linalg.norm(result.centered_positive_energy-result.centered_positive_energy.conj().T)),
                      'full_C1b':'OPEN'},indent=2))

if __name__=='__main__':main()
