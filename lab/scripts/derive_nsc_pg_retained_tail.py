#!/usr/bin/env python3
"""Positive-energy high-energy tail inputs for the other 31 retained massive groups.

Reuses the published group13 MassiveTailPackets/EndpointTailSeries method.
Stores raw positive-energy K,G and Fourier/normalization bounds only.
Does not signed-fold, assemble a covariance, or close full retained C1b.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_pg_tail_packets import MassiveTailPackets
from recursive_horizons.nsc_pg_tail_series import EndpointTailSeries
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-pg-retained-tail.json'
PHYSICAL_SOURCES=(
    'src/recursive_horizons/nsc_pg_tail_packets.py',
    'src/recursive_horizons/nsc_pg_tail_series.py',
    'src/recursive_horizons/nsc_pg_high_energy.py',
    'src/recursive_horizons/nsc_pg_lll_preparation.py',
    'src/recursive_horizons/nsc_lorentzian.py',
    'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
    'src/recursive_horizons/nsc_paired_horizon_preparation.py',
    'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
)
RECORD_SOURCES=PHYSICAL_SOURCES+(
    'scripts/derive_nsc_pg_retained_tail.py',
    'tests/test_nsc_pg_retained_tail.py',
    'docs/nsc-pg-retained-tail.md',
)
INPUTS=(
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-compact-matched-restart.json',
    'results/development/nsc-pg-group13-covariance.json',
    'docs/nsc-pg-tail-moment-audit.json',
)
L0=160.;L1=320.;ORDER=16;GRID=64;POWER=24;NORM_TERMS=6
C_TOL=1e-9;G_TOL=2e-9
MODE_ORDER_REPRESENTATIVES=((1,1),(14,1),(23,1),(32,1))
EXCLUDED=13

def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def sci(value):
    x=float(value)
    if not np.isfinite(x):return None
    return x

def context():
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    restart=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())
    return state,restart['scattering_provenance']['config'],restart['locked_inputs']

def families(state):
    rows=[]
    for ch in state['channels']:
        index=int(ch['index'])
        if index==0 or index==EXCLUDED:continue
        mass=float(ch['compact_mass']);ang=float(ch['angular_eigenvalue'])
        signs=(1,-1) if ang else (1,)
        for sign in signs:rows.append((index,int(sign),mass,sign*ang))
    if len(rows)!=61:raise ValueError(f'expected 61 signed families, got {len(rows)}')
    return rows

def _cache_path(cache,group,sign,lower,order,grid):
    tag=f'g{group:02d}_s{sign:+d}_L{int(lower)}_o{order}_n{grid}.npz'
    return Path(cache)/tag

def evaluate_tail(job):
    group,sign,mass,angular,lower,order,grid,config,cache=job
    path=_cache_path(cache,group,sign,lower,order,grid)
    label=np.array([group,sign,mass,angular,lower,order,grid],float)
    if path.exists():
        with np.load(path,allow_pickle=False) as a:stored={k:a[k].copy() for k in a.files}
        if not np.allclose(stored['label'],label,rtol=0,atol=1e-15):
            raise ValueError('work-cache label differs: '+path.name)
        return stored
    os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
    os.environ['MKL_NUM_THREADS']='1';os.environ['NUMEXPR_NUM_THREADS']='1'
    with threadpool_limits(limits=1):
        try:
            p=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],
                                   config['horizon_offset'],config['scattering_tolerance'])
            packets=MassiveTailPackets(p,mass,angular,order=order,points_per_unit=grid,minimum_energy=lower)
            series=EndpointTailSeries(packets,lower,maximum_power=POWER,normalization_terms=NORM_TERMS)
            K,bK=series.centered_integral();G,bG=series.integral((1.,1.,1.))
            arrays={'label':label,'centered':K,'gram':G,'bounds':np.array([bK,bG],float),
                    'ok':np.array([1],np.uint8),'error':np.frombuffer(b'',np.uint8)}
        except Exception as exc:
            arrays={'label':label,'centered':np.full((8,8),np.nan,complex),'gram':np.full((8,8),np.nan,complex),
                    'bounds':np.array([np.inf,np.inf],float),'ok':np.array([0],np.uint8),
                    'error':np.frombuffer((type(exc).__name__+': '+str(exc)).encode(),np.uint8)}
    path.write_bytes(deterministic_npz_bytes(arrays))
    return arrays

def _pool(workers):
    if workers>2:raise ValueError('at most 2 numerical workers')
    # mpmath precision contexts must never be shared across worker threads.
    return ProcessPoolExecutor(max_workers=workers)

def _run(jobs,workers,label):
    if workers>2:raise ValueError('at most 2 numerical workers')
    out=[]
    print(f'{label}: {len(jobs)} jobs; workers={workers}',flush=True)
    if workers==1:
        return [evaluate_tail(job) for job in jobs]
    with _pool(workers) as pool:
        for n,result in enumerate(pool.map(evaluate_tail,jobs,chunksize=1),1):
            out.append(result)
            if n==len(jobs) or n%8==0:print(f'{label} {n}/{len(jobs)}',flush=True)
    return out

def _key(row):return int(row['label'][0]),int(row['label'][1])

def _group_ok(parts):
    if any(int(p['ok'][0])==0 for p in parts):return False
    c=sum(float(p['bounds'][0]) for p in parts);g=sum(float(p['bounds'][1]) for p in parts)
    return c<=C_TOL and g<=G_TOL

def prepare(workers=2):
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[key]='1'
    state,config,locked=context();rows=families(state)
    source_hashes={p:sha(p) for p in PHYSICAL_SOURCES}
    signature={'sources':source_hashes,'config':config,'method':[L0,L1,ORDER,GRID,POWER,NORM_TERMS]}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-retained-tail-'+tag);cache.mkdir(exist_ok=True)
    jobs=[(g,s,m,a,L0,ORDER,GRID,config,str(cache)) for g,s,m,a in rows]
    first=_run(jobs,workers,f'L={int(L0)} order{ORDER}')
    by_group={}
    for item in first:by_group.setdefault(int(item['label'][0]),[]).append(item)
    need_refine=sorted(g for g,parts in by_group.items() if not _group_ok(parts))
    print('groups requiring L=320: '+json.dumps(need_refine),flush=True)
    refined={}
    if need_refine:
        jobs320=[]
        for g,s,m,a in rows:
            if g in need_refine:jobs320.append((g,s,m,a,L1,ORDER,GRID,config,str(cache)))
        for item in _run(jobs320,workers,f'L={int(L1)} order{ORDER}'):
            refined[_key(item)]=item
    selected=[];diagnostics=[]
    open_families=[]
    group_L={};group_C={};group_G={};group_pass={}
    for g,s,m,a in rows:
        base=next(item for item in first if _key(item)==(g,s))
        if g in need_refine:
            diagnostics.append(base)
            chosen=refined[(g,s)]
        else:
            chosen=base
        selected.append(chosen)
    for g,parts in ((g,[x for x in selected if int(x['label'][0])==g]) for g in sorted({r[0] for r in rows})):
        group_L[g]=float(parts[0]['label'][4])
        if any(float(p['label'][4])!=group_L[g] for p in parts):
            raise ValueError(f'group {g} mixed tail lower bounds')
        group_C[g]=sum(float(p['bounds'][0]) for p in parts)
        group_G[g]=sum(float(p['bounds'][1]) for p in parts)
        group_pass[g]=_group_ok(parts)
        if not group_pass[g]:
            for p in parts:
                open_families.append({'group':int(p['label'][0]),'sign':int(p['label'][1]),
                                      'L':float(p['label'][4]),'ok':bool(int(p['ok'][0])),
                                      'C':sci(p['bounds'][0]),'G':sci(p['bounds'][1]),
                                      'error':bytes(p['error']).decode()})
    checks=[]
    check_jobs=[]
    selected_lookup={_key(item):item for item in selected}
    for g,s in MODE_ORDER_REPRESENTATIVES:
        item=selected_lookup[(g,s)];lower=float(item['label'][4]);mass=float(item['label'][2]);ang=float(item['label'][3])
        check_jobs.append((g,s,mass,ang,lower,12,GRID,config,str(cache)))
    if check_jobs:
        for item in _run(check_jobs,workers,'mode-order 12 vs 16'):
            g,s=_key(item);ref=selected_lookup[(g,s)]
            if int(item['ok'][0]) and int(ref['ok'][0]):
                dC=float(np.linalg.norm(ref['centered']-item['centered'],2))
                dG=float(np.linalg.norm(ref['gram']-item['gram'],2))
            else:
                dC=dG=float('nan')
            checks.append({'group':g,'sign':s,'L':float(ref['label'][4]),'dC':sci(dC),'dG':sci(dG),
                           'order12_ok':bool(int(item['ok'][0])),'order12_C':sci(item['bounds'][0]),
                           'order12_G':sci(item['bounds'][1])})
    arrays={
        'group':np.array([int(x['label'][0]) for x in selected],np.int32),
        'sign':np.array([int(x['label'][1]) for x in selected],np.int32),
        'compact_mass':np.array([float(x['label'][2]) for x in selected],float),
        'signed_angular':np.array([float(x['label'][3]) for x in selected],float),
        'selected_lower':np.array([float(x['label'][4]) for x in selected],float),
        'centered':np.array([x['centered'] for x in selected],complex),
        'gram':np.array([x['gram'] for x in selected],complex),
        'bound_C':np.array([float(x['bounds'][0]) for x in selected],float),
        'bound_G':np.array([float(x['bounds'][1]) for x in selected],float),
        'ok':np.array([int(x['ok'][0]) for x in selected],np.uint8),
        'diag160_group':np.array([int(x['label'][0]) for x in diagnostics],np.int32),
        'diag160_sign':np.array([int(x['label'][1]) for x in diagnostics],np.int32),
        'diag160_centered':np.array([x['centered'] for x in diagnostics],complex) if diagnostics else np.zeros((0,8,8),complex),
        'diag160_gram':np.array([x['gram'] for x in diagnostics],complex) if diagnostics else np.zeros((0,8,8),complex),
        'diag160_bound_C':np.array([float(x['bounds'][0]) for x in diagnostics],float),
        'diag160_bound_G':np.array([float(x['bounds'][1]) for x in diagnostics],float),
        'diag160_ok':np.array([int(x['ok'][0]) for x in diagnostics],np.uint8),
    }
    metadata={'schema':'NSC-PG-RETAINED-TAIL-INPUT-v1','source_hashes':source_hashes,'config':config,
              'locked_inputs':locked,'excluded_group':EXCLUDED,'families':61,
              'signed_fold':False,'full_covariance_assigned':False,
              'method':{'order':ORDER,'points_per_unit':GRID,'maximum_power':POWER,
                        'normalization_terms':NORM_TERMS,'L0':L0,'L1':L1},
              'groups_requiring_320':need_refine,'open_families':open_families,
              'group_L':{str(k):v for k,v in group_L.items()},
              'group_C_sum':{str(k):v for k,v in group_C.items()},
              'group_G_sum':{str(k):v for k,v in group_G.items()},
              'group_pass':{str(k):bool(v) for k,v in group_pass.items()},
              'mode_order_12_vs_16':checks,
              'tail_scope':'positive-energy mode endpoint Fourier series plus explicit Fourier/normalization remainder; no stress subtraction, no signed energy fold, no covariance assembly'}
    arrays['metadata_json']=np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-pg-retained-tail.{digest}.npz';path.write_bytes(raw)
    print(json.dumps({'input':str(path.relative_to(ROOT)),'sha256':digest,'bytes':len(raw),
                      'groups_requiring_320':need_refine,'open_families':len(open_families)},indent=2))
    return path,arrays,metadata

def make_record(path,arrays,metadata):
    open_families=metadata['open_families']
    passed=not open_families
    status=('PASS: tail-only positive-energy high-energy inputs for 31 retained massive groups '
            '(61 signed families); full retained C1b OPEN' if passed else
            'OPEN: some retained-tail families fail the L=320 representation bound; full retained C1b OPEN')
    groups=sorted(int(k) for k in metadata['group_L'])
    record={'schema':'NSC-PG-RETAINED-TAIL-v1',
            'status':status,
            'accountable_author':'Douglas Ek',
            'source_hashes':{p:sha(p) for p in RECORD_SOURCES},
            'input_hashes':{p:sha(p) for p in INPUTS},
            'payload':{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                       'bytes':path.stat().st_size},
            'locked_inputs':metadata['locked_inputs'],
            'precision_audit':{'evidence':'docs/nsc-pg-tail-moment-audit.json','tables':29,'maximum_numerical_change':0.0,'audit_decimal_precision':50,'future_execution':'separate processes or explicit serial mode; no thread fallback'},
            'domain':{'surface':'same common PG Cauchy surface on the fixed transmitted Dirac background',
                      'energy':'positive-energy tail only; negative energy is not assigned by a partner fold',
                      'excluded_group':EXCLUDED,
                      'families':61,
                      'angular0_plus_only':[23],
                      'signed_fold':False,
                      'full_covariance_assigned':False},
            'method':{**metadata['method'],
                      'packets':'MassiveTailPackets(order=16, points_per_unit=64, minimum_energy=L)',
                      'series':'EndpointTailSeries(lower=L, maximum_power=24, normalization_terms=6)',
                      'centered':'centered_integral()',
                      'gram':'integral((1,1,1))',
                      'L_is_numerical_split':True},
            'tolerances':{'C_bound_sum':C_TOL,'Gram_bound_sum':G_TOL},
            'groups_requiring_320':metadata['groups_requiring_320'],
            'open_families':open_families,
            'coverage':{'groups':groups,'family_count':int(len(arrays['group'])),
                        'angular0_plus_only':[23],'excluded':[EXCLUDED],
                        'selected_lower':{str(g):int(metadata['group_L'][str(g)]) for g in groups},
                        'C_bound_sum':{str(g):sci(metadata['group_C_sum'][str(g)]) for g in groups},
                        'G_bound_sum':{str(g):sci(metadata['group_G_sum'][str(g)]) for g in groups},
                        'max_C_bound_sum':sci(max(metadata['group_C_sum'].values())),
                        'max_G_bound_sum':sci(max(metadata['group_G_sum'].values())),
                        'all_pass':passed},
            'mode_order_12_vs_16':metadata['mode_order_12_vs_16'],
            'error_account':{'representation':'Fourier/normalization remainder of the positive-energy tail',
                             'mode_approximation':'order 12 vs 16 at selected L; not a uniform physical remainder theorem',
                             'seed_covariance_used_as_input':False,
                             'LLL_covariance_used_as_substitute':False},
            'gate':{'tail_inputs':'PASS' if passed else 'OPEN',
                    'full_C1b':'OPEN',
                    'group13_covariance':'PASS',
                    'signed_fold':False,
                    'full_covariance_assigned':False,
                    'physical_stress_computed':False,
                    'metric_timestep':False,
                    'PDF_bump':False,
                    'A_q_Omega_zeta_Vfull_changed':False},
            'comparison':{'fields':'all JSON fields and NPZ array identities',
                          'routine':'authenticate the tail artifact; no old generator rerun'}}
    text=json.dumps(record,indent=2,sort_keys=True)+'\n'
    return record,text

def load_artifact(path):
    raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    if path.name!=f'nsc-pg-retained-tail.{digest}.npz':raise ValueError('tail artifact digest differs')
    with np.load(path,allow_pickle=False) as a:arrays={k:a[k].copy() for k in a.files}
    metadata=json.loads(arrays['metadata_json'].tobytes())
    return arrays,metadata,digest

def verify_arrays(arrays,metadata,state):
    rows=families(state)
    if list(map(int,arrays['group']))!=[r[0] for r in rows]:raise ValueError('group coverage differs')
    if list(map(int,arrays['sign']))!=[r[1] for r in rows]:raise ValueError('sign coverage differs')
    if EXCLUDED in set(map(int,arrays['group'])):raise ValueError('group13 must remain excluded')
    for g,s,m,a in rows:
        if g==23 and s!=1:raise ValueError('angular0 group23 must be plus only')
    for g,parts in ((g,[(i,int(arrays['sign'][i])) for i in range(len(arrays['group'])) if int(arrays['group'][i])==g])
                    for g in sorted(set(map(int,arrays['group'])))):
        lowers={float(arrays['selected_lower'][i]) for i,_ in parts}
        if len(lowers)!=1:raise ValueError(f'group {g} mixed selected L')
    for i,(g,s,m,a) in enumerate(rows):
        if arrays['compact_mass'][i]!=m or arrays['signed_angular'][i]!=a:
            raise ValueError('tail physical label differs')
        if int(arrays['ok'][i]):
            if not all(np.isfinite(arrays[key][i]).all() for key in ('centered','gram','bound_C','bound_G')):
                raise ValueError('selected tail input is nonfinite')
            if np.linalg.norm(arrays['gram'][i]-arrays['gram'][i].conj().T)>3e-11:
                raise ValueError('tail Gram is not Hermitian')
    for g in sorted(set(map(int,arrays['group']))):
        idx=np.flatnonzero(arrays['group']==g)
        C=float(np.sum(arrays['bound_C'][idx]));G=float(np.sum(arrays['bound_G'][idx]))
        if not np.isclose(C,metadata['group_C_sum'][str(g)],rtol=3e-13,atol=0) or not np.isclose(G,metadata['group_G_sum'][str(g)],rtol=3e-13,atol=0):
            raise ValueError('stored group error budget differs from arrays')
        passed=bool(np.all(arrays['ok'][idx])) and C<=C_TOL and G<=G_TOL
        if passed!=metadata['group_pass'][str(g)]:raise ValueError('group gate differs from actual bounds')
    if metadata['signed_fold'] or metadata['full_covariance_assigned']:
        raise ValueError('tail artifact must not assign a folded covariance')

def authenticate(path,arrays,metadata):
    state,_,locked=context()
    if locked!=metadata['locked_inputs']:raise ValueError('locked A/Omega/zeta/V_full inputs changed')
    for name,expected in metadata['source_hashes'].items():
        if sha(name)!=expected:raise ValueError('physical input dependency changed: '+name)
    verify_arrays(arrays,metadata,state)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--prepare-inputs',action='store_true')
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--input')
    parser.add_argument('--record',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.workers>2:raise ValueError('at most 2 numerical workers')
    if args.check:
        stored=json.loads((ROOT/OUTPUT).read_text())
        path=ROOT/stored['payload']['path'];arrays,metadata,digest=load_artifact(path)
        authenticate(path,arrays,metadata)
        record,_=make_record(path,arrays,metadata)
        if stored['payload']['sha256']!=digest:raise ValueError('stored payload hash differs')
        if stored!=json.loads(json.dumps(record,sort_keys=True)):
            raise ValueError('stored tail record differs from authenticated recombination')
        print(json.dumps({'status':stored['status'],'gate':stored['gate'],
                          'groups_requiring_320':stored['groups_requiring_320']},indent=2));return
    if args.prepare_inputs:
        prepare(args.workers);return
    if args.record:
        if not args.input:raise ValueError('explicit authenticated tail artifact required')
        path=ROOT/args.input;arrays,metadata,_=load_artifact(path)
        authenticate(path,arrays,metadata)
        record,text=make_record(path,arrays,metadata)
        (ROOT/OUTPUT).write_text(text)
        print(json.dumps({'status':record['status'],'bytes':len(text.encode()),
                          'payload':record['payload'],'gate':record['gate']},indent=2));return
    path,arrays,metadata=prepare(args.workers)
    record,text=make_record(path,arrays,metadata)
    (ROOT/OUTPUT).write_text(text)
    print(json.dumps({'status':record['status'],'bytes':len(text.encode()),
                      'payload':record['payload'],'groups_requiring_320':record['groups_requiring_320'],
                      'gate':record['gate']},indent=2))

if __name__=='__main__':main()
