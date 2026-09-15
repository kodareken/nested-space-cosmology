#!/usr/bin/env python3
"""New massive PG mode fields, with current/Green-jump and outer-error records.

No historical generator is invoked. Existing tables supply labels and
comparison controls; the spatial fields come from the owned Dirac equation.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_pg_massive_modes import MassivePGModeResolution,physical_sewing
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

CODE=('scripts/derive_nsc_pg_massive_modes.py','src/recursive_horizons/nsc_massive_jost_modes.py',
      'src/recursive_horizons/nsc_pg_massive_modes.py')
DEPENDENCIES=('src/recursive_horizons/nsc_paired_horizon_preparation.py',
              'src/recursive_horizons/nsc_unruh_state.py','src/recursive_horizons/nsc_gauge_source.py',
              'src/recursive_horizons/nsc_pg_lll_preparation.py','src/recursive_horizons/nsc_lorentzian.py',
              'src/recursive_horizons/nsc_common_time_bulk_split.py','src/recursive_horizons/nsc_transmitting_dirac_domain.py',
              'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py')
INNER=np.array([-1.,0.,1.]);OUTER=np.array([2.1,3.,4.])

def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def inputs():
    matched=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    for p in (matched['scattering_payload'],matched['seed_generations_payload'],state['payload']):
        if sha(p['path'])!=p['sha256']:raise ValueError('authenticated input differs: '+p['path'])
    # Only declared frequency/weight arrays are read from the old seed here.
    with np.load(ROOT/state['payload']['path'],allow_pickle=False) as a:
        frequencies=a['frequency'].copy();weights=a['quadrature_weight'].copy()
    labels=[]
    for ch in state['channels'][1:]:
        for sign in ((1,-1) if ch['angular_eigenvalue'] else (1,)):
            for index in range(ch['sample_offset'],ch['sample_offset']+ch['sample_count']):
                labels.append([ch['index'],sign,index,ch['compact_mass'],sign*ch['angular_eigenvalue'],frequencies[index],weights[index]])
    return np.array(labels),matched,state

def worker(job):
    i,row,config,cache,refine=job;path=Path(cache)/f'{i:04d}.npz'
    if path.exists():
        with np.load(path,allow_pickle=False) as a:
            if not np.array_equal(a['label'],row):raise ValueError('work-cache labels differ')
        return i
    p=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],config['horizon_offset'],config['scattering_tolerance'])
    owner=MassivePGModeResolution(p)
    _,_,_,m,lam,E,_=row
    jost=solve_jost(p.background,E,m,lam,radial_collar=config['horizon_offset'])
    modes=owner.section(E,m,lam,INNER,OUTER,jost=jost)
    residuals={**modes.residuals,**modes.spectral_jump_residual(),**jost.residuals}
    inner,_,future,projection=physical_sewing(E,m,jost.reflection,jost.transmission)
    scattering=np.zeros((3,3),complex);scattering[:2]=inner
    if projection[2,2]:scattering[2]=jost.transmission_phase*future[2]
    residuals['complex_scattering_isometry']=float(np.linalg.norm(scattering.conj().T@scattering-projection))
    residuals['complex_t_current']=float(abs(abs(jost.complex_transmission)**2-jost.transmission))
    if refine:
        farther=solve_jost(p.background,E,m,lam,radial_collar=config['horizon_offset'],outer_radius=2*jost.outer_radius)
        refined=owner.section(E,m,lam,INNER,OUTER,jost=farther)
        residuals['outer_radius_doubling_reflection']=float(abs(jost.reflection-farther.reflection))
        residuals['outer_radius_doubling_PG_fields']=float(max(np.linalg.norm(modes.interior-refined.interior),np.linalg.norm(modes.exterior-refined.exterior)))
    eig=np.linalg.eigvalsh(modes.source_covariance)
    residuals['source_CAR']=max(0.,-float(eig.min()),float(eig.max())-1.)
    frame=owner.frame(0.)[0]
    modal=np.linalg.solve(frame,modes.interior[1])
    arrays={'label':row,'interior':modes.interior,'exterior':modes.exterior,
            'interior_fundamental':modes.interior_fundamental,'exterior_ingoing':modes.exterior_ingoing,
            'source_covariance':modes.source_covariance,'source_projector':projection,
            'scattering':scattering,'reflection':np.array(jost.reflection),'transmission':np.array(jost.transmission),
            'complex_transmission':np.array(jost.complex_transmission),
            'seed_control_covariance':modal@modes.source_covariance@modal.conj().T,
            'residuals_json':np.frombuffer(json.dumps(residuals,sort_keys=True).encode(),np.uint8)}
    path.write_bytes(deterministic_npz_bytes(arrays))
    return i

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--limit',type=int);parser.add_argument('--cache-only',action='store_true')
    args=parser.parse_args();labels,matched,state=inputs();config=matched['scattering_provenance']['config']
    signature=hashlib.sha256(json.dumps({p:sha(p) for p in (*CODE,*DEPENDENCIES)},sort_keys=True).encode()).hexdigest()
    cache=Path(tempfile.gettempdir())/f'nsc-pg-massive-modes-{signature[:16]}';cache.mkdir(exist_ok=True)
    # One outer-radius refinement per class at a below/above-threshold node.
    refined=set()
    for ch in (1,13,14,22,23,32):
        candidates=np.flatnonzero((labels[:,0]==ch)&(labels[:,1]==1))
        for fraction in (.13,.65):refined.add(int(candidates[int(fraction*len(candidates))]))
    work=list(range(len(labels)))[:args.limit]
    print(f'New PG mode fields: {len(work)} labels; cache {cache}',flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for n,i in enumerate(pool.map(worker,((i,labels[i],config,str(cache),i in refined) for i in work),chunksize=1),1):
            if n%64==0 or n==len(work):print(f'PG mode fields {n}/{len(work)}',flush=True)
    if args.cache_only or args.limit:
        print('Finite worker cache only; no full-state or complete-grid claim.',flush=True);return
    results=[];data=[]
    for i in range(len(labels)):
        with np.load(cache/f'{i:04d}.npz',allow_pickle=False) as a:
            data.append({k:a[k].copy() for k in a.files if k!='residuals_json'})
            results.append(json.loads(a['residuals_json'].tobytes()))
    arrays={k:np.array([d[k] for d in data]) for k in data[0]}
    arrays['rho_interior']=INNER;arrays['rho_exterior']=OUTER
    # Controls are loaded only after physical mode/source evaluation completes.
    with np.load(ROOT/matched['seed_generations_payload']['path'],allow_pickle=False) as a:
        old_labels=a['labels'].copy();old_matched=a['covariance_matched_delta_q'].copy()
    lookup={(int(row[0]),int(row[1]),int(row[2])):j for j,row in enumerate(old_labels)}
    comparison=[]
    for ch in state['channels'][1:]:
        mask=labels[:,0]==ch['index'];indices=np.flatnonzero(mask)
        diffs=[]
        if ch['compact_level']:
            for i in indices:
                old_i=lookup[tuple(labels[i,:3].astype(int))]
                diffs.append(float(np.linalg.norm(arrays['seed_control_covariance'][i]-old_matched[old_i])))
        maxima={k:max(results[i].get(k,0.) for i in indices) for k in results[0]}
        comparison.append({'channel_index':ch['index'],'positive_energy_signed_angular_nodes':len(indices),
                           'compact_mass':ch['compact_mass'],'residual_maxima':maxima,
                           'matched_seed_control_difference':max(diffs) if diffs else None})
    metadata={'schema':'NSC-PG-MASSIVE-MODE-FIELDS-v1','source_hashes':{p:sha(p) for p in (*CODE,*DEPENDENCIES)},
              'label_columns':['channel','angular_sign','seed_control_index','mass','signed_angular','energy','weight'],
              'state_input':'C_H plus inherited incoming occupation; seed covariance control only',
              'spectral_measure':'dE/(2*pi); samples are not an identity-CAR finite discretization',
              'normalization':'complex Jost amplitudes and matched_delta_q inner modes',
              'prior_seed_generation':matched['start_conventions']['matched_delta_q']['version']}
    arrays['metadata_json']=np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)
    binary=deterministic_npz_bytes(arrays);digest=hashlib.sha256(binary).hexdigest()
    payload=ROOT/f'results/development/artifacts/nsc-pg-massive-mode-fields.{digest}.npz';payload.write_bytes(binary)
    all_keys=set().union(*(r.keys() for r in results))
    summary={'source_hashes':metadata['source_hashes'],'locked_inputs':matched['locked_inputs'],
             'payload':{'path':str(payload.relative_to(ROOT)),'sha256':digest,'bytes':len(binary)},
             'nodes':len(labels),'groups':comparison,
             'maximum_residuals':{k:max(r.get(k,0.) for r in results) for k in sorted(all_keys)},
             'refinements':[{'label_index':i,'label':labels[i].tolist(),'residuals':results[i]} for i in sorted(refined)]}
    out=cache/'evaluation.json';out.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'evaluation':str(out),'payload':summary['payload'],'maximum_residuals':summary['maximum_residuals']},indent=2))

if __name__=='__main__':main()
