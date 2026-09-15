#!/usr/bin/env python3
"""Evaluate the same subgap state with bounded retarded packet kernels."""
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

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_pg_threshold_projection as previous
from recursive_horizons.nsc_complex_horizon_modes import finite_pg_horizon_basis,complex_reflection
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_pg_packet_modes import project_horizon_bases,SIGNED_PACKET_MAP
from recursive_horizons.nsc_pg_retarded_packets import retarded_packets
from recursive_horizons.nsc_pg_retarded_threshold import retarded_threshold_integral
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes


def evaluate(job):
    z,config,mass,angular,old_cache,cache=job;key=previous.node_key(z);path=Path(cache)/(key+'.npz')
    if path.exists():return
    p=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],config['horizon_offset'],config['scattering_tolerance'])
    basis=finite_pg_horizon_basis(z,mass,angular,p)
    old_path=Path(old_cache)/(key+'.npz')
    if old_path.exists():
        with np.load(old_path,allow_pickle=False) as a:
            arrays={k:a[k].copy() for k in a.files if k!='residuals_json'}
        if arrays['frequency'].item()!=z:raise ValueError('imported mode frequency differs')
    else:
        F,_=project_horizon_bases(z,mass,angular,basis.interior,basis.exterior,rtol=2e-13,atol=2e-15)
        Fc=F.copy() if not z.imag else project_horizon_bases(z.conjugate(),mass,angular,basis.interior_conjugate,basis.exterior_conjugate,rtol=2e-13,atol=2e-15)[0]
        arrays={'frequency':np.array(z),'projection':F,'projection_conjugate':Fc}
        if z.imag:
            arrays['reflection']=np.array(complex_reflection(z,mass,angular,p.background,radial_collar=config['horizon_offset']).reflection)
    residuals={'basis':basis.residuals}
    if z.imag:
        G=retarded_packets(z,mass,angular,p,basis=basis,projected_basis=arrays['projection'])
        arrays['packet_resolvent']=G.response;residuals['packet_resolvent']=G.residuals
    arrays['residuals_json']=np.frombuffer(json.dumps(residuals,sort_keys=True).encode(),np.uint8)
    path.write_bytes(deterministic_npz_bytes(arrays))


def main():
    p=argparse.ArgumentParser();p.add_argument('--channel',type=int,default=13);p.add_argument('--sign',type=int,choices=(-1,1),default=1)
    p.add_argument('--points',type=int,default=24);p.add_argument('--height-fraction',type=float,default=3.);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--adaptive',action='store_true');p.add_argument('--epsabs',type=float,default=1e-10)
    args=p.parse_args()
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    config=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    channel=next(c for c in state['channels'] if c['index']==args.channel)
    mass=float(channel['compact_mass']);angular=args.sign*float(channel['angular_eigenvalue'])
    if mass<=1 or (not angular and args.sign!=1) or args.height_fraction<=2:raise ValueError('existing compact channel and pole-free contour required')
    old_signature={'source_hashes':{q:previous.sha(q) for q in previous.DEPENDENCIES},
                   'node_function_sha256':hashlib.sha256(inspect.getsource(previous.evaluate_node).encode()).hexdigest(),
                   'config':config,'channel':channel,'angular_sign':args.sign}
    old_tag=hashlib.sha256(json.dumps(old_signature,sort_keys=True).encode()).hexdigest()[:16]
    old_cache=Path(tempfile.gettempdir())/('nsc-threshold-contour-'+old_tag)
    signature={**old_signature,'retarded_sources':{q:previous.sha(q) for q in
               ('src/recursive_horizons/nsc_pg_retarded_packets.py','src/recursive_horizons/nsc_pg_retarded_threshold.py')},
               'retarded_node_function_sha256':hashlib.sha256(inspect.getsource(evaluate).encode()).hexdigest()}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-retarded-threshold-'+tag);cache.mkdir(exist_ok=True)
    if args.adaptive:
        from recursive_horizons.nsc_pg_adaptive_threshold import adaptive_threshold_integral
        memo={}
        def node(z):
            z=complex(z)
            if z not in memo:
                evaluate((z,config,mass,angular,str(old_cache),str(cache)))
                with np.load(cache/(previous.node_key(z)+'.npz'),allow_pickle=False) as a:
                    memo[z]={k:a[k].copy() for k in a.files if k!='residuals_json'}
                if len(memo)%32==0:print(f'Adaptive channel {args.channel} sign {args.sign}: {len(memo)} nodes',flush=True)
            return memo[z]
        print(f'Adaptive bounded subgap: channel {args.channel}, sign {args.sign}; cache {cache}',flush=True)
        result=adaptive_threshold_integral(1.,mass,config['surface_gravity'],node,
                                           height=config['surface_gravity']/args.height_fraction,epsabs=args.epsabs)
        data={'gram_positive':result.gram_positive_energy,'centered_positive':result.centered_positive_energy,
              'analytic':result.analytic_integrals,'smooth':result.smooth_integral}
        if not angular:
            data['gram_signed']=(result.gram_positive_energy+SIGNED_PACKET_MAP@result.gram_positive_energy.conj()@SIGNED_PACKET_MAP)/(2*np.pi)
            data['centered_signed']=(result.centered_positive_energy-SIGNED_PACKET_MAP@result.centered_positive_energy.conj()@SIGNED_PACKET_MAP)/(2*np.pi)
        metadata={**signature,'integration_source_sha256':previous.sha('src/recursive_horizons/nsc_pg_adaptive_threshold.py'),
                  'adaptive':True,'epsabs':args.epsabs,'height_fraction':args.height_fraction,'interval':[1.,mass],
                  'diagnostics':result.diagnostics,'evaluated_nodes':len(memo),
                  'scope':'subgap panel only; adaptive errors are numerical estimates, not a global tail bound'}
        data['metadata_json']=np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)
        path=cache/f'panel-adaptive-h{args.height_fraction:g}-tol{args.epsabs:g}.npz'
        path.write_bytes(deterministic_npz_bytes(data))
        print(json.dumps({'panel':str(path),'diagnostics':result.diagnostics,
                          'Gram_positive_eigenvalues':np.linalg.eigvalsh(result.gram_positive_energy).tolist(),
                          'full_C1b':'OPEN'},indent=2))
        return
    nodes=list(dict.fromkeys(previous.nodes_for(1.,mass,config['surface_gravity'],args.points,args.height_fraction)))
    print(f'Bounded subgap response: channel {args.channel}, sign {args.sign}, {len(nodes)} nodes; cache {cache}',flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i,_ in enumerate(pool.map(evaluate,((z,config,mass,angular,str(old_cache),str(cache)) for z in nodes),chunksize=1),1):
            if i%16==0 or i==len(nodes):print(f'Bounded nodes {i}/{len(nodes)}',flush=True)
    values={}
    for z in nodes:
        with np.load(cache/(previous.node_key(z)+'.npz'),allow_pickle=False) as a:
            values[z]={k:a[k].copy() for k in a.files if k!='residuals_json'}
    result=retarded_threshold_integral(1.,mass,config['surface_gravity'],lambda z:values[z],points=args.points,
                                      height=config['surface_gravity']/args.height_fraction)
    data={'gram_positive':result.gram_positive_energy,'centered_positive':result.centered_positive_energy,
          'analytic':result.analytic_integrals,'smooth':result.smooth_integral}
    if not angular:
        data['gram_signed']=(result.gram_positive_energy+SIGNED_PACKET_MAP@result.gram_positive_energy.conj()@SIGNED_PACKET_MAP)/(2*np.pi)
        data['centered_signed']=(result.centered_positive_energy-SIGNED_PACKET_MAP@result.centered_positive_energy.conj()@SIGNED_PACKET_MAP)/(2*np.pi)
    metadata={**signature,'points':args.points,'height_fraction':args.height_fraction,'interval':[1.,mass],
              'scope':'bounded representation of same subgap panel; full energy tail remains unassigned'}
    data['metadata_json']=np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)
    path=cache/f'panel-p{args.points}-h{args.height_fraction:g}.npz';path.write_bytes(deterministic_npz_bytes(data))
    print(json.dumps({'panel':str(path),'Gram_positive_trace':float(np.trace(data['gram_positive']).real),
                      'Gram_positive_eigenvalues':np.linalg.eigvalsh(data['gram_positive']).tolist(),
                      'Hermitian_residual':float(np.linalg.norm(data['centered_positive']-data['centered_positive'].conj().T)),
                      'full_C1b':'OPEN'},indent=2))

if __name__=='__main__':main()
