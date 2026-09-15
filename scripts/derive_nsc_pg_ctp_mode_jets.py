#!/usr/bin/env python3
"""Bind physical PG source modes to the transmitting history insertion."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_ctp_mode_jets import continue_mode_columns,mode_metric_vertices,PGSourceModeHistoryDerivative
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_dirac_domain import S3,MODE_TO_CURRENT
from recursive_horizons.nsc_pg_lll_preparation import PGLLLPreparation,LLLFlowAtlas
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_lorentzian import geometry

OUTPUT='results/development/nsc-pg-ctp-mode-jets.json'
SOURCES=('src/recursive_horizons/nsc_pg_ctp_mode_jets.py','scripts/derive_nsc_pg_ctp_mode_jets.py',
         'tests/test_nsc_pg_ctp_mode_jets.py','docs/nsc-pg-ctp-mode-jets.md')
INPUTS=('results/development/nsc-pg-retained-covariance.json','results/development/nsc-pg-massive-mode-resolution.json',
        'results/development/nsc-pg-lll-preparation.json','results/development/nsc-transmitting-ctp-resolvent.json',
        'results/development/nsc-mode-resolved-cauchy-state.json','results/development/nsc-compact-matched-restart.json',
        'src/recursive_horizons/nsc_transmitting_resolvent.py','src/recursive_horizons/nsc_transmitting_ctp_variation.py',
        'src/recursive_horizons/nsc_common_time_bulk_split.py','scripts/derive_nsc_transmitting_boundary_binding.py')

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def rec(p):return json.loads((ROOT/p).read_text())
def load(p):
    with np.load(p,allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}

def input_families():
    mode=rec(INPUTS[1]);mp=mode['payload']
    if sha(mp['path'])!=mp['sha256']:raise ValueError('mode artifact changed')
    a=load(ROOT/mp['path']);channels=rec(INPUTS[4])['channels'];config=rec(INPUTS[5])['scattering_provenance']['config']
    families=[]
    for ch in channels[1:]:
        g=ch['index'];sgns=(1,-1) if ch['angular_eigenvalue'] else (1,)
        for sign in sgns:
            i=np.flatnonzero((a['label'][:,0]==g)&(a['label'][:,1]==sign))
            chosen=i[[int(np.argmin(abs(a['label'][i,5]-target))) for target in (.2,.6)]]
            E=a['label'][chosen,5];w=a['label'][chosen,6];P=a['source_projector'][chosen];C=a['source_covariance'][chosen]
            other_sign=-sign if ch['angular_eigenvalue'] else sign
            other=[]
            for energy in E:
                j=np.flatnonzero((a['label'][:,0]==g)&(a['label'][:,1]==other_sign)&(a['label'][:,5]==energy))
                if len(j)!=1:raise ValueError('physical opposite-angular mode is absent')
                other.append(int(j[0]))
            Fplus=a['interior'][chosen,1]
            Fminus=np.einsum('ij,njk->nik',S3,a['interior'][other,1].conj())
            families.append({'group':g,'sign':sign,'mass':ch['compact_mass'],'angular':sign*ch['angular_eigenvalue'],
              'energies':np.r_[-E[::-1],E],'weights':np.r_[w[::-1],w],
              'at_zero':np.concatenate((Fminus[::-1],Fplus)),
              'source':np.concatenate((P[::-1]-C[::-1].conj(),C)),'projector':np.concatenate((P[::-1],P)),
              'indices':np.r_[np.array(other)[::-1],chosen]})
    # LLL uses its own exact characteristic map, never assigned to masses.
    state=rec(INPUTS[4]);a0=load(ROOT/state['payload']['path']);ch=channels[0]
    sl=slice(ch['sample_offset'],ch['sample_offset']+ch['sample_count']);oldE=a0['frequency'][sl]
    selected=[int(np.argmin(abs(oldE-x))) for x in (.2,.6)];E=oldE[selected]
    # Frequency weights are solely a finite differential-control rule.
    # Reuse the declared angular frequency grid weights from the state payload.
    weight_key=next((k for k in ('weight','frequency_weight','quadrature_weight') if k in a0),None)
    w=a0[weight_key][sl][selected] if weight_key else np.ones(2)
    lll=PGLLLPreparation(config['horizon_rho'],config['surface_gravity'],config['omega'])
    beta=float(geometry(0.)[0]);F=np.zeros((2,3),complex)
    F[:,1]=MODE_TO_CURRENT[:,0]/np.sqrt(beta-1);F[:,2]=MODE_TO_CURRENT[:,1]/np.sqrt(beta+1)
    signed=np.r_[-E[::-1],E]
    families.insert(0,{'group':0,'sign':1,'mass':0.,'angular':0.,'energies':signed,'weights':np.r_[w[::-1],w],
          'at_zero':np.repeat(F[None,:,:],len(signed),axis=0),
          'source':np.array([lll.spectral_covariance(e) for e in signed]),
          'projector':np.repeat(np.eye(3)[None,:,:],len(signed),axis=0),
          'indices':np.full(4,-1,dtype=int)})
    return families

def one(job):
    family,cache=job;g,s=family['group'],family['sign'];out=Path(cache)/f'group-{g}-sign-{s}.npz'
    if out.exists():return str(out)
    arrays={k:np.asarray(v) for k,v in family.items()};meta={}
    for points in (48,72):
        fields=continue_mode_columns(family['energies'],family['mass'],family['angular'],family['at_zero'],points=points)
        M,D=mode_metric_vertices(fields,family['mass'],family['angular'])
        arrays[f'vertex_{points}']=M;arrays[f'direct_{points}']=D
        meta[f'seam_{points}']=fields['seam_residual']
        if points==72:
            for k,v in fields.items():arrays['fields/'+k]=np.asarray(v)
    arrays['metadata_json']=np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)
    out.write_bytes(deterministic_npz_bytes(arrays));return str(out)

def prepare(workers):
    inputs={p:sha(p) for p in INPUTS};source={p:sha(p) for p in SOURCES[:2]}
    signature={'source':source,'input':inputs};tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-pg-ctp-mode-jets-'+tag);cache.mkdir(exist_ok=True)
    families=input_families();print(f'{len(families)} physical signed families; {cache}',flush=True)
    arrays={}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for n,path in enumerate(pool.map(one,((f,str(cache)) for f in families)),1):
            d=load(path);name=f"{int(d['group'].item())}_{int(d['sign'].item())}"
            for k,v in d.items():arrays[name+'/'+k]=v
            if n%8==0:print(f'{n}/{len(families)}',flush=True)
    arrays['metadata_json']=np.frombuffer(json.dumps({'signature':signature,'families':len(families),'physical_history_selected':False},sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest()
    out=ROOT/f'results/development/artifacts/nsc-pg-ctp-mode-jets.{digest}.npz';out.write_bytes(raw);print(out,flush=True)
    return out

def evaluate(a):
    names=sorted({k.split('/')[0] for k in a if '/' in k});rows={};maxima={}
    for name in names:
        d={k[len(name)+1:]:v for k,v in a.items() if k.startswith(name+'/')}
        M=d['vertex_72'];D=d['direct_72'];owner=PGSourceModeHistoryDerivative(d['energies'],M,d['source'],d['projector'])
        controls=owner.source_fiber_ctp_control(d['weights'])
        row={'group':int(d['group'].item()),'sign':int(d['sign'].item()),'energies':d['energies'].tolist(),
          'weak_vs_direct':float(np.max(abs(M-D))),'spatial_refinement':float(np.max(abs(M-d['vertex_48']))),
          'energy_exchange_hermiticity':float(np.max(abs(M-M.swapaxes(-1,-2).conj()))),
          'source_CAR_lower_violation':max(0.,-float(np.linalg.eigvalsh(d['source']).min())),
          'source_CAR_upper_violation':max(0.,-float(np.linalg.eigvalsh(d['projector']-d['source']).min())),
          'CTP_tangent':max(x['tangent_residual'] for x in controls),
          'CTP_trace':max(x['source_trace_residual'] for x in controls),
          'CTP_linear_solve':max(x['linear_solve_residual'] for x in controls),
          'seam':float(d['fields/seam_residual'].item()),
          'source_rank':controls[0]['open_source_rank'],
          'positive_frequency_transfer_norm':float(np.linalg.norm(M[:,6:9,9:12])),
          'source_differential_control':[float(x['derivative'].real) for x in controls]}
        rows[name]=row
        for k in ('weak_vs_direct','spatial_refinement','energy_exchange_hermiticity','source_CAR_lower_violation',
                  'source_CAR_upper_violation','CTP_tangent','CTP_trace','CTP_linear_solve','seam'):
            maxima[k]=max(maxima.get(k,0.),row[k])
    if {r['group'] for r in rows.values()}!=set(range(33)):raise ValueError('not all retained groups evaluated')
    return rows,maxima

def make_record(path):
    a=load(path);metadata=json.loads(a['metadata_json'].tobytes())
    for group in metadata['signature'].values():
        for p,digest in group.items():
            if sha(p)!=digest:raise ValueError('physical input owner changed: '+p)
    rows,maxima=evaluate(a);tol=3e-8;algebra=3e-11
    failures={k:v for k,v in maxima.items() if v>(tol if k in ('weak_vs_direct','spatial_refinement') else algebra)}
    return {'schema':'NSC-PG-CTP-MODE-JETS-v1','accountable_author':'Douglas Ek',
      'status':'PASS: physical PG mode-space history insertion; full KS endpoint jets OPEN' if not failures else 'OPEN',
      'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
      'payload':{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size},
      'locked_inputs':rec(INPUTS[0])['locked_inputs'],
      'domain':{'surface':'fixed common-PG source-mode representation, rho=0 transmitting seam',
        'metric_fields':['log_N','beta','log_q_PG','log_r'],'spatial_direction':'existing smooth compact profile on (-1,1)',
        'spectral_measure':'dE/(2*pi); no discrete sample is a complete Cauchy space',
        'state':'same C_H plus inherited incoming occupation; opposite-energy covariance on open source fibers',
        'new_particle_sectors':False,'packet_covariance_as_C0':False,'given_B':False,
        'seam_account':'full continuous fields; no sharp-projection delta distribution relabeled as link or stress'},
      'construction':{'vertex':'M_A(Eo,Ei)=integral Phi_Eo dagger deltaH_A Phi_Ei d rho',
        'history_derivative':'deltaU/dg_A(tau)=-i exp[-i Eo(tf-tau)] M_A exp[-i Ei(tau-ti)]',
        'source_control':'finite Nyström differential through existing ctp_first_variation; not an absolute stress or full action trace',
        'kernel_scope':'first variation at the owned static PG background; times remain arguments, no history selected',
        'whole_bulk_retained':True,'ordinary_matrix_link_assigned':False},
      'per_family':rows,'maxima':maxima,'tolerances':{'spatial_and_differential':tol,'algebra':algebra},'failures':failures,
      'jet_inventory':{'PG_source_mode_metric_insertion':'computed','PG_time_ordered_derivative_kernel':'defined from computed insertion',
        'CTP_source_fiber_contraction':'computed control','raw_KS_endpoint_pullback':None,'transmitting_EndpointBranchJets':None},
      'gate':{'B1_full':'OPEN: raw KS endpoint/Cauchy-slice pullback and history-dependent domain must be supplied',
        'Gamma_rest':'OPEN','physical_two_sided_Weyl_mismatch':None,'Weyl_diagnostic':93.54264532195464,
        'stationarity':'OPEN','history_selected':False,'metric_timestep':False,'stress_computed':False,'Z3':'OUT OF SCOPE'},
      'comparison':{'fields':'all','float_atol':3e-12,'float_rtol':3e-10,'exact':'source/input hashes and domain labels'},
      'reproducer':'python3 scripts/derive_nsc_pg_ctp_mode_jets.py --check'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--workers',type=int,default=2)
    p.add_argument('--input');p.add_argument('--check',action='store_true');args=p.parse_args()
    path=prepare(args.workers) if args.prepare else (ROOT/rec(OUTPUT)['payload']['path'] if args.check else Path(args.input).resolve())
    record=make_record(path)
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(rec(OUTPUT),record)
    else:(ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':record['status'],'maxima':record['maxima'],'failures':record['failures']},indent=2))

if __name__=='__main__':main()
