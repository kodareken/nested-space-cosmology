#!/usr/bin/env python3
"""Apply the actual PG-to-KS coordinate map to computed transmitting modes."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_ks_metric_pullback import pulled_mode_vertices,coordinate_residuals,reference_chart,pg_log_jacobian,fixed_surface_normal
from recursive_horizons.nsc_pg_ctp_mode_jets import PGSourceModeHistoryDerivative
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-pg-ks-metric-pullback.json'
INPUT='results/development/nsc-pg-ctp-mode-jets.json'
SOURCES=('src/recursive_horizons/nsc_pg_ks_metric_pullback.py','scripts/derive_nsc_pg_ks_metric_pullback.py',
         'tests/test_nsc_pg_ks_metric_pullback.py','docs/nsc-pg-ks-metric-pullback.md')
DEPENDENCIES=(INPUT,'src/recursive_horizons/nsc_pg_ctp_mode_jets.py','src/recursive_horizons/nsc_lorentzian.py',
 'src/recursive_horizons/nsc_common_time_bulk_split.py','src/recursive_horizons/nsc_transmitting_dirac_domain.py',
 'src/recursive_horizons/nsc_transmitting_resolvent.py','docs/nsc-adm-neck-source-map.md',
 'src/recursive_horizons/nsc_transmitting_boundary_history.py','scripts/derive_nsc_transmitting_boundary_binding.py')
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def rec(p):return json.loads((ROOT/p).read_text())
def load(p):
    with np.load(p,allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}
def one(job):
    name,fields,mass,angular,cache=job;path=Path(cache)/(name+'.npz')
    if not path.exists():
        M,D=pulled_mode_vertices(fields,mass,angular)
        path.write_bytes(deterministic_npz_bytes({'raw_KS_vertex':M,'direct_vertex':D}))
    return name,str(path)

def prepare(workers):
    r=rec(INPUT)
    if sha(r['payload']['path'])!=r['payload']['sha256']:raise ValueError('PG source-mode input changed')
    a=load(ROOT/r['payload']['path']);signature={p:sha(p) for p in (*SOURCES[:2],*DEPENDENCIES)}
    tag=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:16]
    cache=Path(tempfile.gettempdir())/('nsc-pg-ks-pullback-'+tag);cache.mkdir(exist_ok=True)
    names=sorted({k.split('/')[0] for k in a if '/' in k});jobs=[]
    for name in names:
        d={k[len(name)+1:]:v for k,v in a.items() if k.startswith(name+'/')}
        fields={k[len('fields/'):]:v for k,v in d.items() if k.startswith('fields/')}
        jobs.append((name,fields,float(d['mass'].item()),float(d['angular'].item()),str(cache)))
    arrays={}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for count,(name,path) in enumerate(pool.map(one,jobs),1):
            for k,v in load(path).items():arrays[name+'/'+k]=v
            if count%16==0:print(f'{count}/{len(names)} existing mode families',flush=True)
    arrays['metadata_json']=np.frombuffer(json.dumps({'signature':signature,'input_payload':r['payload'],'new_mode_integration':False},sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest()
    out=ROOT/f'results/development/artifacts/nsc-pg-ks-metric-pullback.{digest}.npz';out.write_bytes(raw);print(out,flush=True)
    return out

def evaluate(a,mode):
    rows={};maxima={}
    for name in sorted({k.split('/')[0] for k in a if '/' in k}):
        d={k[len(name)+1:]:v for k,v in mode.items() if k.startswith(name+'/')}
        M=a[name+'/raw_KS_vertex'];D=a[name+'/direct_vertex'];g,s=map(int,name.split('_'))
        state=PGSourceModeHistoryDerivative(d['energies'],M,d['source'],d['projector'])
        controls=state.source_fiber_ctp_control(d['weights'])
        _,_,a0,_,_=reference_chart(0.);J0=pg_log_jacobian(0.,1.,0.,a0,1.)
        incorrect_constant=np.einsum('AB,Aij->Bij',J0,d['vertex_72'])
        row={'group':g,'sign':s,'weak_vs_direct':float(np.max(abs(M-D))),
          'Hermiticity':float(np.max(abs(M-M.swapaxes(-1,-2).conj()))),
          'CTP_tangent':max(x['tangent_residual'] for x in controls),
          'CTP_trace':max(x['source_trace_residual'] for x in controls),
          'constant_neck_matrix_error_control':float(np.linalg.norm(M-incorrect_constant)),
          'source_differential_control':[float(x['derivative'].real) for x in controls]}
        rows[name]=row
        for k in ('weak_vs_direct','Hermiticity','CTP_tangent','CTP_trace'):
            maxima[k]=max(maxima.get(k,0.),row[k])
    chart={str(x):coordinate_residuals(x) for x in (-.8,0.,.8)}
    _,_,a0,_,_=reference_chart(0.);normal=fixed_surface_normal(0.,1.,0.,a0,1.)
    normal={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in normal.items()}
    return rows,maxima,chart,normal

def make_record(path):
    a=load(path);metadata=json.loads(a['metadata_json'].tobytes())
    for p,digest in metadata['signature'].items():
        if sha(p)!=digest:raise ValueError('pullback dependency changed: '+p)
    mode=load(ROOT/metadata['input_payload']['path'])
    rows,maxima,chart,normal=evaluate(a,mode)
    failures={k:v for k,v in maxima.items() if v>(3e-8 if k=='weak_vs_direct' else 3e-11)}
    for point,values in chart.items():
        for k,v in values.items():
            if v>(3e-8 if k=='jacobian_vs_finite_difference' else 3e-11):failures[point+'/'+k]=v
    _,_,a0,_,J=reference_chart(0.);jac=pg_log_jacobian(0.,1.,0.,a0,1.)
    return {'schema':'NSC-PG-KS-METRIC-PULLBACK-v1','accountable_author':'Douglas Ek',
      'status':'PASS: raw KS metric directions on the reference chart; two-endpoint history pullback OPEN' if not failures else 'OPEN',
      'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in DEPENDENCIES},
      'payload':{'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size},
      'locked_inputs':rec(INPUT)['locked_inputs'],
      'domain':{'chart':'dT=-d rho/a_ref, dz=d tau+beta_ref*d rho/a_ref^2; chart held fixed during metric variation',
        'reference':'a_ref=sqrt(beta_ref^2-1) on the already owned trapped geometry',
        'raw_KS_fields':['N','beta','q_ADM','r'],'q_ADM':'a_parallel, not q_PG',
        'variation_profile':'existing s(rho) applied to raw KS metric coefficients through the reference chart',
        'time_dependent_directions':'may use the conditional PG-time insertion; this is not a homogeneous KS two-endpoint interpolation',
        'normal':'fixed rho=0 normal metric variation; not a Cauchy-state isometry',
        'moving_surface':False,'new_metric_solution':False},
      'seed_event':{'coordinate_jacobian':J.tolist(),'PG_log_by_KS_raw_jacobian':jac.tolist(),'determinant':float(np.linalg.det(jac)),'fixed_normal':normal},
      'chart_controls':chart,'per_family':rows,'maxima':maxima,'failures':failures,
      'tolerances':{'coordinate_and_CTP_algebra':3e-11,'independent_finite_difference':3e-8},
      'jet_inventory':{'PG_source_mode_history_kernel':'imported computed','raw_KS_metric_profile_kernel':'computed with rho-dependent Jacobian',
        'fixed_surface_normal_derivative':'computed at the neck','KS_two_endpoint_shape_and_embedding':None,'full_EndpointBranchJets':None},
      'gate':{'B1_full':'OPEN: actual two-endpoint history/Cauchy embedding remains to be supplied',
        'Gamma_rest':'OPEN','physical_two_sided_Weyl_mismatch':None,'Weyl_diagnostic':93.54264532195464,
        'stationarity':'OPEN','history_selected':False,'metric_timestep':False,'stress_computed':False,'Z3':'OUT OF SCOPE'},
      'comparison':{'fields':'all','float_atol':3e-13,'float_rtol':3e-13,'exact':'source/input hashes and domain labels'},
      'reproducer':'python3 scripts/derive_nsc_pg_ks_metric_pullback.py --check'}

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
