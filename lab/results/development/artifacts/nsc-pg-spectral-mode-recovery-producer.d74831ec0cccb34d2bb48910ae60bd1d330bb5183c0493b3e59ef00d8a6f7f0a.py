#!/usr/bin/env python3
"""Recover missing real spectral fields from authenticated packet columns."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_mode_recovery import stationary_interior_observation, StationaryInteriorObservation
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-pg-spectral-mode-recovery.json'
SOURCES=('src/recursive_horizons/nsc_pg_mode_recovery.py','tests/test_nsc_pg_mode_recovery.py',
         'scripts/derive_nsc_pg_spectral_mode_recovery.py','docs/nsc-pg-spectral-mode-recovery.md')
INPUTS=('results/development/nsc-pg-retained-covariance.json','results/development/nsc-pg-group13-covariance.json',
        'results/development/nsc-pg-massive-mode-resolution.json','results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-compact-matched-restart.json','src/recursive_horizons/nsc_pg_retarded_packets.py',
        'src/recursive_horizons/nsc_pg_lll_preparation.py','src/recursive_horizons/nsc_lorentzian.py',
        'src/recursive_horizons/nsc_paired_horizon_preparation.py','src/recursive_horizons/nsc_transmitting_dirac_domain.py')

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
    with np.load(p,allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}
def authenticated(p):
    if sha(p['path'])!=p['sha256']:raise ValueError('payload changed: '+p['path'])
    return load(ROOT/p['path'])

def jobs():
    retained=read(INPUTS[0]); a=authenticated(retained['payload'])
    meta=json.loads(a['metadata_json'].tobytes());g13=read(INPUTS[1]);b=authenticated(g13['payload'])
    c=authenticated(meta['group14_input']);channels=read(INPUTS[3])['channels']
    cfg=read(INPUTS[4])['scattering_provenance']['config'];out=[]
    def append(name, group, sign, E, w, F, C=None, P=None):
        E=np.asarray(E).ravel();w=np.asarray(w).ravel();ch=channels[group];m=ch['compact_mass']
        expected_C=np.array([source_covariance(e,cfg['surface_gravity'],cfg['omega'],m) for e in E])
        expected_P=np.array([np.diag([1.,1.,float(e>m)]) for e in E])
        if C is not None and np.max(abs(C-expected_C))>3e-13:raise ValueError('source law changed '+name)
        if P is not None and np.max(abs(P-expected_P))>3e-13:raise ValueError('open source fiber changed '+name)
        out.append((name,group,sign,m,sign*ch['angular_eigenvalue'],E,w,F,expected_C,expected_P))
    prefixes=sorted({k.rsplit('/',1)[0] for k in a if k.endswith('/projection') and k.startswith(('low/','low_ref/','mid/','mid_ref/'))})
    for prefix in prefixes:
        low=prefix.startswith('low');group,sign=map(int,prefix.split('/')[1].split('_'))
        get=lambda key:a[prefix+'/'+key]
        append(prefix,group,sign,get('energy' if low else 'energies'),get('weight' if low else 'weights'),
               get('projection'),get('source') if low else None,get('projector') if low else None)
    append('group13/low',13,1,b['low_energy'],b['low_weight'],b['low_projection'],b['low_source'],b['low_projector'])
    append('group13/mid',13,1,b['fast_energy_24'],b['fast_weight_24'],b['fast_projection_24'])
    for sign in (1,-1):
        for prefix in (f'low16_{sign}',f'low32_{sign}',f'mid24_{sign}'):
            low=prefix.startswith('low');get=lambda key:c[prefix+'/'+key]
            append('group14/'+prefix,14,sign,get('energy' if low else 'energies'),get('weight' if low else 'weights'),
                   get('projection'),get('source') if low else None,get('projector') if low else None)
    return out,[retained['payload'],g13['payload'],meta['group14_input']]

def one(job):
    name,group,sign,m,l,E,w,F,C,P=job
    owner=stationary_interior_observation(E,m,l)
    result=owner.recover(F,P,require_all=False)
    return name,{'group':np.array(group),'sign':np.array(sign),'mass':np.array(m),'angular':np.array(l),
                 'energies':E,'weights':w,'source':C,'projector':P,'projection':F,
                 'observation':owner.observation,'fundamental_at_zero':owner.fundamental_at_zero,
                 'fundamental_at_minus_one':owner.fundamental_at_minus_one,'current_residual':owner.current_residual,
                 'function_evaluations':np.array(owner.function_evaluations),
                 'mode_at_one':result.mode_at_one,'accepted':result.accepted}

def prepare(workers):
    work,provenance=jobs(); arrays={};signature={p:sha(p) for p in (*SOURCES,*INPUTS)}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for n,(name,row) in enumerate(pool.map(one,work),1):
            arrays.update({name+'/'+k:v for k,v in row.items()})
            if n%10==0 or n==len(work):print(f'{n}/{len(work)} real panels reconstructed; no horizon/scattering rerun',flush=True)
    if signature!={p:sha(p) for p in (*SOURCES,*INPUTS)}:raise ValueError('dependency changed during preparation')
    arrays['metadata_json']=np.frombuffer(json.dumps({'signature':signature,'input_payloads':provenance},sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-pg-spectral-mode-recovery.{digest}.npz';path.write_bytes(raw)
    return path

def make_record(path):
    a=load(path);meta=json.loads(a['metadata_json'].tobytes())
    for p,s in meta['signature'].items():
        if sha(p)!=s:raise ValueError('recovery dependency changed '+p)
    for p in meta['input_payloads']:
        if sha(p['path'])!=p['sha256']:raise ValueError('source payload changed')
    panels={};total=accepted=0;maxima={};groups=set();failure=[]
    names=sorted(k[:-len('/mode_at_one')] for k in a if k.endswith('/mode_at_one'))
    for name in names:
        get=lambda k:a[name+'/'+k]
        owner=StationaryInteriorObservation(get('energies'),float(get('mass')),float(get('angular')),get('observation'),
                 get('fundamental_at_zero'),get('fundamental_at_minus_one'),get('current_residual'),int(get('function_evaluations')))
        result=owner.recover(get('projection'),get('projector'),require_all=False)
        if not np.array_equal(result.accepted,get('accepted')):raise ValueError('acceptance mask changed '+name)
        if not np.allclose(result.mode_at_one,get('mode_at_one'),rtol=3e-11,atol=3e-12,equal_nan=True):raise ValueError('recovered fields changed '+name)
        ok=result.accepted;E=get('energies');group=int(get('group'));groups.add(group)
        rowmax={k:max((d[k] for d,b in zip(result.diagnostics,ok) if b and d[k] is not None),default=None)
                for k in ('input_error_amplification','projection_absolute_residual','current_coisometry_residual','source_projector_residual','transport_current_residual')}
        for k,v in rowmax.items():
            if v is not None:maxima[k]=max(maxima.get(k,0.),v)
        bad=[{'index':i,'energy':float(E[i]),'status':d['status'],'sigma_min':d['singular_values'][1],
              'input_error_amplification':d['input_error_amplification']} for i,d in enumerate(result.diagnostics) if not ok[i]]
        panels[name]={'group':group,'sign':int(get('sign')),'rows':len(E),'accepted':int(ok.sum()),
                      'energy_range':[float(E.min()),float(E.max())],
                      'accepted_energy_range':[float(E[ok].min()),float(E[ok].max())] if ok.any() else None,
                      'rejected':bad,'maxima_on_accepted':rowmax,
                      'status':'PASS' if ok.all() else 'OPEN: unresolved rows explicitly excluded'}
        total+=len(E);accepted+=int(ok.sum())
    if groups!=set(range(1,33)):failure.append('massive group coverage incomplete')
    return {'schema':'NSC-PG-SPECTRAL-MODE-RECOVERY-v1','accountable_author':'Douglas Ek',
      'status':'PASS: resolved real-frequency mode recovery; full spectral source remains OPEN' if not failure else 'OPEN',
      'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
      'payload':{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size},
      'input_payloads':meta['input_payloads'],'locked_inputs':read(INPUTS[0])['locked_inputs'],
      'panels':panels,'counts':{'panels':len(panels),'groups':len(groups),'rows':total,'accepted_rows':accepted,'unresolved_rows':total-accepted},
      'maxima_on_accepted':maxima,'failures':failure,
      'tolerances':{'projection':3e-11,'field':3e-8,'current':3e-8,'declared_input_accuracy':3e-11},
      'domain':{'input':'authenticated frequency-resolved physical packet columns, never integrated covariance or seed',
        'operator':'same PG Dirac radial ODE on [-1,1]; identity fundamental convention at rho=1',
        'inverse':'resolved six-by-two observation matrix; no singular-value floor',
        'source':'unchanged C_H plus inherited incoming occupation, exact open fiber',
        'panel_accounting':'base and refinement panels remain separate; no overlapping panels added here',
        'negative_energy':'already owned opposite-angular sigma3 K map; positive real source rows stored',
        'scope':'conditional input fields for the common KS source, not a finite closed packet dynamics',
        'horizon_scattering_rerun':False,'state_changed':False,'history_selected':False},
      'gate':{'resolved_mode_recovery':'PASS' if not failure else 'OPEN','all_real_rows':'PASS' if accepted==total else 'OPEN',
        'full_spectral_source':'OPEN: missing subgap/tail fields and common nonlinear subtraction/convergence',
        'stress':None,'nulls':None,'updated_constraints':None,'metric_timestep':False,'Z3':'OUT OF SCOPE',
        'V_c':None,'extended_stationarity':'OPEN','PDF_bumped':False,'Weyl_time_node_diagnostic':93.54264532195464,
        'homogeneous_nonexistence':'preserved; not rerun'},
      'comparison':{'fields':'all','float_atol':3e-12,'float_rtol':3e-10,'exact':'hashes, labels, masks, coverage and scope'},
      'reproducer':'python3 scripts/derive_nsc_pg_spectral_mode_recovery.py --check'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--check',action='store_true');p.add_argument('--workers',type=int,default=3);p.add_argument('--input');args=p.parse_args()
    path=prepare(args.workers) if args.prepare else (ROOT/read(OUTPUT)['payload']['path'] if args.check else Path(args.input).resolve())
    result=make_record(path)
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT),result)
    else:(ROOT/OUTPUT).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','counts','maxima_on_accepted','failures')},indent=2),flush=True)
    if result['failures']:raise SystemExit(1)

if __name__=='__main__':main()
