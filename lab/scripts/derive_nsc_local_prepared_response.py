#!/usr/bin/env python3
"""Four bounded source-fixed nonzero-history runs on the resolved local grid.

The new gap is actual nonzero-history/PDE matter and tangent accuracy. Reuse
the7601 initial columns and all existing sources. Center/plus/minus at60
steps check the full derivative; center at120 steps isolates temporal error.
At most four solves and300 CPU seconds; --check uses saved restrictions only.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_compatible_history_geometry import CompatibleIncomingMetric,CompatibleRadiusDirection
from recursive_horizons.nsc_evolved_incoming_state import AmplitudeOnlyMetric,FixedSourcePreparation
from recursive_horizons.nsc_exact_phase_prepared_state import prepare_exact_phase_incoming
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter,surface_geometry_response,compatible_history_slots
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator
from derive_nsc_evolved_incoming_constraints import coefficients_from_records

OUTPUT='results/development/nsc-local-prepared-response.json'
ARTIFACT='results/development/artifacts/nsc-local-prepared-response.npz'
OWNED=('scripts/derive_nsc_local_prepared_response.py','docs/nsc-local-prepared-response.md')
INPUTS=('results/development/nsc-local-reference-accuracy.json',
        'results/development/nsc-incoming-source-update-v5.json',
        'results/development/nsc-incoming-surface-coefficients.json',
        'results/development/nsc-incoming-surface-regular-branch.json',
        'results/development/nsc-incoming-fixed-transfer.json',
        'src/recursive_horizons/nsc_exact_phase_prepared_state.py',
        'src/recursive_horizons/nsc_evolved_incoming_state.py',
        'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
        'src/recursive_horizons/nsc_local_incoming_family.py',
        'src/recursive_horizons/nsc_compatible_history_geometry.py',
        'src/recursive_horizons/nsc_transmitting_history_jets.py',
        'src/recursive_horizons/nsc_transmitting_history_modes.py',
        'scripts/derive_nsc_evolved_incoming_constraints.py')
ALPHA=.001
STEP=.0001
CAP=300.
CASES=(('center',ALPHA,60,True),('plus',ALPHA+STEP,60,False),
       ('minus',ALPHA-STEP,60,False),('time_fine',ALPHA,120,True))


class BudgetReached(RuntimeError):pass
def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}
def norm(x):return float(np.max(abs(x),initial=0.))


def family(alpha):
    center=chart_coordinates(1.)[1]+.15
    w=LocalAxialFunction((0.,1.),center)
    zero=LocalAxialFunction((0.,),center)
    return CompatibleIncomingMetric((alpha,),(CompatibleRadiusDirection(w,zero,.007,.03),))


def load_inputs():
    records=[json.loads((ROOT/p).read_text()) for p in INPUTS[:5]]
    p=records[0]['payload']
    if digest(p['path'])!=p['sha256']:raise ValueError('refined initial column artifact changed')
    with np.load(ROOT/p['path'],allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    channel=records[4]['actual_source_bindings']['channel']
    if channel['compact_mass']!=float(arrays['source/mass'].item()) or channel['angular_eigenvalue']!=float(arrays['source/angular'].item()):
        raise ValueError('source/channel operator changed')
    return records,arrays,channel


def evaluate(arrays,name,alpha,coeff,baseline,parameters):
    def get(key):return arrays[name+'/'+key]
    times=get('times');z=times+chart_coordinates(1.)[1]
    F,Fz,dF,dFz=(get(k) for k in ('F','Fz','dF','dFz'))
    weights=arrays['source/column_weights'];C=arrays['source/covariance'];E=arrays['source/energies_repeated']
    current=source_column_matter(F*weights,Fz*weights,C,dF*weights,dFz*weights,**parameters)
    ref=arrays['reference_F0'][None]*np.exp(-1j*times[:,None,None]*E)
    refz=-1j*E*ref
    zero=np.zeros((0,*ref.shape),complex)
    reference=source_column_matter(ref*weights,refz*weights,C,zero,zero,**parameters)
    slots,ds=compatible_history_slots(family(alpha),z,len(dF))
    geom=surface_geometry_response(slots,ds,coeff)
    value=baseline+geom['action_gradient_change']+current['action_gradient']-reference['action_gradient']
    tangent=geom['action_gradient_tangent']+current['action_gradient_tangent']
    return {'residual':value,'tangent':tangent,'matter':current['action_gradient'],
            'matter_tangent':current['action_gradient_tangent'],
            'matter_correction':current['action_gradient']-reference['action_gradient']}


def analyze(arrays,meta):
    rows={}
    for name,alpha,steps,_ in CASES:
        if name not in meta['completed']:continue
        got=evaluate(arrays,name,alpha,meta['coefficients'],np.asarray(meta['baseline']),meta['parameters'])
        mask=np.arange(steps+1)>=2*steps//3
        rows[name]={'sampled_local_constraint_approximant_max':np.max(abs(got['residual'][mask]),axis=0).tolist(),
                    'sampled_local_matter_correction_max':np.max(abs(got['matter_correction'][mask]),axis=0).tolist(),
                    'phase_residual':meta['cases'][name]['phase_residual'],
                    'flux_residual':meta['cases'][name]['flux_residual']}
    checks={}
    if all(n in rows for n in ('center','plus','minus')):
        center=evaluate(arrays,'center',ALPHA,meta['coefficients'],np.asarray(meta['baseline']),meta['parameters'])
        plus=evaluate(arrays,'plus',ALPHA+STEP,meta['coefficients'],np.asarray(meta['baseline']),meta['parameters'])
        minus=evaluate(arrays,'minus',ALPHA-STEP,meta['coefficients'],np.asarray(meta['baseline']),meta['parameters'])
        fd=(plus['residual']-minus['residual'])/(2*STEP)
        checks['complete_local_constraint_derivative']=norm((fd-center['tangent'][0])[40:])
        checks['local_matter_derivative']=norm(((plus['matter']-minus['matter'])/(2*STEP)-center['matter_tangent'][0])[40:])
        checks['dropped_deltaC_effect']=norm(center['matter_tangent'][0,40:])
        checks['same_fixed_preparation']=len({meta['cases'][n]['preparation'] for n in ('center','plus','minus')})==1
        checks['derivative_passed']=checks['complete_local_constraint_derivative']<=3e-8 and checks['same_fixed_preparation']
    if 'time_fine' in rows and 'center' in rows:
        coarse=evaluate(arrays,'center',ALPHA,meta['coefficients'],np.asarray(meta['baseline']),meta['parameters'])
        fine=evaluate(arrays,'time_fine',ALPHA,meta['coefficients'],np.asarray(meta['baseline']),meta['parameters'])
        checks['local_time_change']=np.max(abs(fine['matter'][::2][40:]-coarse['matter'][40:]),axis=0).tolist()
        checks['local_tangent_time_change']=np.max(abs(fine['matter_tangent'][0,::2][40:]-coarse['matter_tangent'][0,40:]),axis=0).tolist()
        checks['time_indicator_target']=1e-11
        checks['time_indicator_passed']=max(checks['local_time_change'])<=1e-11
    return {'per_case':rows,'checks':checks}


def record(arrays,meta,payload):
    measurements=analyze(arrays,meta)
    return {'schema':'NSC-LOCAL-PREPARED-RESPONSE-v1','accountable_author':'Douglas Ek',
            'status':'OPEN: source-fixed nonzero-history control; full-source local gate uncertified',
            'control':{'amplitude':ALPHA,'finite_difference_step':STEP,'runs_max':4,
                       'profile':'w=alpha*Chebyshev_T1*flat axial cutoff, U=0; normal .007/.03',
                       'family':'14_1','points':7601,'local_interval':'S(1)+[.12,.18]',
                       'derivative_tolerance':3e-8,'temporal_indicator_tolerance':1e-11},
            'measurements':measurements,'runtime':meta['runtime'],
            'scope':{'physical_local_gate':'OPEN','full_source_error':None,'field_and_derivative_error':None,
                     'between_node_error':None,'all_channel_coverage':False,'source_changed':False,
                     'stress_drift_subtracted':False,'frozen_C0_imposed':False,'metric_timestep':False},
            'source_hashes':{p:meta['signature'][p] for p in OWNED},
            'input_hashes':{p:meta['signature'][p] for p in INPUTS},
            'reused_artifacts':meta['reused_artifacts'],'payload':payload,
            'reproducer':'python3 scripts/derive_nsc_local_prepared_response.py --check'}


def write(arrays,meta):
    raw=deterministic_npz_bytes({**arrays,'metadata_json':np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)})
    (ROOT/ARTIFACT).write_bytes(raw)
    out=record(arrays,meta,{'path':ARTIFACT,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    (ROOT/OUTPUT).write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return out


def run():
    if (ROOT/OUTPUT).exists() or (ROOT/ARTIFACT).exists():raise FileExistsError('existing control; use --check')
    sig=signatures();records,a,channel=load_inputs();coeff=coefficients_from_records(records[2],records[3])
    source=FixedSourcePreparation(a['source/covariance'],a['source/column_weights'],np.repeat(a['source/energies'],3))
    meta={'signature':sig,'completed':[],'attempted':[],'cases':{},'coefficients':{k:np.asarray(v).tolist() for k,v in coeff.items()},
          'baseline':records[1]['baseline']['action_gradient_approximant'],
          'parameters':{'mass':channel['compact_mass'],'angular':channel['angular_eigenvalue'],
                        'axial_scale':coeff['a'],'radius':coeff['r'],'multiplicity':channel['copy_count']*channel['degeneracy']/2},
          'reused_artifacts':[records[0]['payload']]}
    arrays={k:v for k,v in a.items() if k.startswith('source/')}
    arrays['source/energies_repeated']=source.energies
    arrays['reference_F0']=a['reference_columns'][0]
    cpu=time.process_time();wall=time.monotonic();exhausted=False
    def stop(*_):raise BudgetReached('nonzero local prepared300 CPU-second cap')
    def checkpoint():
        meta['runtime']={'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,
                         'CPU_cap':CAP,'attempted':list(meta['attempted']),'completed':list(meta['completed']),'cap_reached':exhausted}
        if signatures()!=sig:raise ValueError('control owner/input changed')
        return write(arrays,meta)
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    try:
        owner=FourthOrderModePropagator(a['x'],channel['compact_mass'],channel['angular_eigenvalue'])
        for name,alpha,steps,tangent in CASES:
            begin=time.process_time();provider=family(alpha)
            if not tangent:provider=AmplitudeOnlyMetric(provider)
            times=np.linspace(0.,.18,steps+1);meta['attempted'].append(name)
            result=prepare_exact_phase_incoming(owner,times,provider,a['initial_fields'],source)
            arrays.update({name+'/times':times,name+'/F':result.state.columns,name+'/Fz':result.axial_columns,
                           name+'/dF':result.state.column_tangents,name+'/dFz':result.axial_tangents})
            meta['cases'][name]={'CPU_seconds':time.process_time()-begin,
                                'preparation':result.fixed_preparation_digest,
                                'history':result.sampled_history_fingerprint,
                                'phase_residual':result.diagnostics['exact_harmonic_phase_residual'],
                                'flux_residual':result.diagnostics['norm_flux_residual']}
            meta['completed'].append(name);del result
            checkpoint();print(name+' complete CPU '+str(meta['cases'][name]['CPU_seconds']),flush=True)
    except BudgetReached:exhausted=True
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    return checkpoint()


def check():
    saved=json.loads((ROOT/OUTPUT).read_text());load_inputs()
    if digest(ARTIFACT)!=saved['payload']['sha256']:raise ValueError('control payload changed')
    with np.load(ROOT/ARTIFACT,allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    if meta['signature']!=signatures():raise ValueError('control source/input changed')
    if record(arrays,meta,saved['payload'])!=saved:raise ValueError('array-only replay differs')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    out=run() if p.parse_args().run else check()
    print(json.dumps({k:out[k] for k in ('status','measurements','runtime')},indent=2))
