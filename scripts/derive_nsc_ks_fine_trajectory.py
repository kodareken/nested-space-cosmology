#!/usr/bin/env python3
"""One fine source-fixed KS trajectory, stored in bounded deterministic chunks."""
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
import derive_nsc_ks_trajectory_control as T
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid,usual_axial_support,physical_incoming_interval,trigonometric_polynomial
from recursive_horizons.nsc_ks_streaming_trajectory import stream_ks_trajectory,TrajectoryChunkWriter,TrajectoryChunkReader
from recursive_horizons.nsc_ks_profile_identity import profile_description,profile_identity
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT=ROOT/'results/development/nsc-ks-fine-trajectory.json'
DIRECTORY=ROOT/'results/development/artifacts/nsc-ks-fine-trajectory'
FINAL=ROOT/'results/development/artifacts/nsc-ks-fine-prepared.npz'
OWNED=('scripts/derive_nsc_ks_fine_trajectory.py','src/recursive_horizons/nsc_ks_streaming_trajectory.py',
       'src/recursive_horizons/nsc_ks_profile_identity.py','docs/nsc-ks-fine-trajectory.md')
INPUTS=tuple(dict.fromkeys((*T.OWNED,*T.INPUTS,T.OUTPUT)))
N,LENGTH,CPU_CAP=2048,205/512,180.
OPTIONS={'rtol':5e-14,'atol':5e-19,'max_step':1/8192}


def digest(path):
    path=Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}


def descriptor(path):
    path=Path(path)
    return {'path':str(path.relative_to(ROOT)),'sha256':digest(path),'bytes':path.stat().st_size}


def save_json(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
    os.replace(temporary,path)


def run():
    if OUTPUT.exists() or FINAL.exists() or DIRECTORY.exists():
        raise FileExistsError('fine trajectory has already been attempted; inspect it instead of repeating')
    sig=signatures();old=json.loads((ROOT/T.OUTPUT).read_text())
    if digest(old['payload']['path'])!=old['payload']['sha256']:raise ValueError('upstream source artifact changed')
    with np.load(ROOT/old['payload']['path'],allow_pickle=False) as f:
        arrays={k:f[k] for k in f.files if k.startswith('source/') or k in ('initial_canonical_columns','target_z')}
        prior_meta=json.loads(f['metadata_json'].tobytes())
    source=FixedSourcePreparation(arrays['source/covariance'],arrays['source/column_weights'],arrays['source/energies'])
    family=T.C.R.C.P.family(T.C.R.C.P.ALPHA)
    grid=computational_z_grid(N,LENGTH)
    if not np.all(np.diff(grid)==LENGTH/N) or grid[0]+LENGTH/2!=family.directions[0].w.center:
        raise ArithmeticError('fine numerical grid is not exactly uniform/centered')
    writer=TrajectoryChunkWriter(DIRECTORY)
    binding={'signature':sig,'parameters':prior_meta['parameters'],'rho_up':prior_meta['rho_up'],
             'options':OPTIONS,'nodes':N,'period_length':LENGTH,'period_origin':float(grid[0]),
             'physical_interval':list(physical_incoming_interval()),'axial_support':list(usual_axial_support()),
             'history_description':profile_description(family),'history_identity':profile_identity(family),
             'axial_profile_identity':profile_identity(family,include_normal_window=False),
             'expected_preparation_digest':prior_meta['preparation'],'reused_artifact':old['payload'],
             'maximum_field_runs':1,'CPU_cap':CPU_CAP,'complete':False}
    save_json(DIRECTORY/'run-binding.json',binding)
    arrays['computational_z']=grid
    cpu,wall=time.process_time(),time.monotonic()
    completed=False;error=None;prepared=None;summary=None
    def stop(*_):raise RuntimeError('fine trajectory180 CPU-second cap')
    def emit(segment):
        writer(segment)
        if writer.buffered_segments==0:
            count=writer._manifest['accepted_steps']
            if count%32==0:print(f'{count} accepted steps saved; rho={segment.rho_end:.9f}',flush=True)
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CPU_CAP)
    try:
        prepared,summary=stream_ks_trajectory(source,arrays['initial_canonical_columns'],family,grid,arrays['target_z'],
            prior_meta['parameters']['mass'],prior_meta['parameters']['angular'],prior_meta['rho_up'],
            axial_support=usual_axial_support(),tangents='zero',on_segment=emit,**OPTIONS)
        if prepared.fixed_preparation_digest!=binding['expected_preparation_digest']:
            raise ValueError('fixed upstream preparation changed')
        prepared.require_history(family)
        completed=True
    except Exception as exc:
        error=f'{type(exc).__name__}: {exc}'
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
        writer.close(complete=completed)
    binding.update(complete=completed,error=error,runtime={'CPU_seconds':time.process_time()-cpu,
        'wall_seconds':time.monotonic()-wall,'CPU_cap':CPU_CAP,'field_runs':1})
    if signatures()!=sig:raise ValueError('fine trajectory implementation/input changed')
    if completed:
        for key in T.FIELDS:arrays['prepared/'+key]=getattr(prepared,key)
        binding.update(preparation=prepared.fixed_preparation_digest,history_fingerprint=prepared.binding.fingerprint,
            accepted_steps=summary.accepted_steps,diagnostics=dict(prepared.diagnostics))
        arrays['metadata_json']=np.frombuffer(json.dumps(binding,sort_keys=True).encode(),np.uint8)
        FINAL.write_bytes(deterministic_npz_bytes(arrays))
    save_json(DIRECTORY/'run-binding.json',binding)
    reader=TrajectoryChunkReader(DIRECTORY,require_complete=completed)
    chunks=[{'path':str((DIRECTORY/e['filename']).relative_to(ROOT)),
             'sha256':e['sha256'],'bytes':e['bytes'],'segment_count':e['segment_count']}
             for e in reader.manifest['chunks']]
    result={'schema':'NSC-KS-FINE-TRAJECTORY-v1','accountable_author':'Douglas Ek',
        'status':('OPEN: fine trajectory captured; continuous error certificate pending' if completed else
                  'OPEN: fine trajectory incomplete; inspect stored chunks'),
        'complete_capture':completed,'error':error,'runtime':binding['runtime'],
        'grid':{'points':N,'period_length_exact':'205/512','binary64_uniform':True,
                'physical_interval':'S(1)+[.12,.18]','physical_profiles_changed':False},
        'options':OPTIONS,'history_identity':binding['history_identity'],
        'axial_profile_identity':binding['axial_profile_identity'],
        'source_preparation_digest':binding['expected_preparation_digest'],
        'trajectory_manifest':descriptor(DIRECTORY/'manifest.json'),
        'run_binding':descriptor(DIRECTORY/'run-binding.json'),'trajectory_chunks':chunks,
        'payload':descriptor(FINAL) if completed else None,
        'source_hashes':{p:sig[p] for p in OWNED},'input_hashes':{p:sig[p] for p in INPUTS},
        'reused_artifact':old['payload'],'scope':{'physical_local_gate':'OPEN','value_error_capture_only':True,
            'full_retarded_derivative_owner_unchanged':True,'source_error_bound':None,
            'continuous_field_error_bound':None,'metric_timestep':False,'stress_drift_subtracted':False},
        'reproducer':'python3 scripts/derive_nsc_ks_fine_trajectory.py --check'}
    save_json(OUTPUT,result)
    return result


def check():
    saved=json.loads(OUTPUT.read_text())
    if not saved['complete_capture']:raise ValueError('fine trajectory is incomplete')
    for path,expected in {**saved['source_hashes'],**saved['input_hashes']}.items():
        if digest(path)!=expected:raise ValueError('fine trajectory source/input binding changed')
    for item in (saved['trajectory_manifest'],saved['run_binding'],saved['payload'],saved['reused_artifact']):
        if digest(item['path'])!=item['sha256']:raise ValueError('fine trajectory payload changed')
    reader=TrajectoryChunkReader(DIRECTORY)
    first,last=None,None;steps=0
    for chunk in reader:
        if first is None:first=np.array(chunk.state_nodes[0],copy=True)
        last=np.array(chunk.state_nodes[-1],copy=True);steps+=len(chunk)
    with np.load(FINAL,allow_pickle=False) as f:a={k:f[k] for k in f.files}
    meta=json.loads(a['metadata_json'].tobytes());ns=len(a['source/energies']);size=2*ns
    if reader.manifest['rho_start']!=meta['rho_up'] or reader.manifest['rho_end']!=1.:
        raise ValueError('fine trajectory does not cover the prepared slab')
    if steps!=meta['accepted_steps'] or len(first)!=size*(1+N):raise ValueError('fine trajectory shape/count mismatch')
    if not np.array_equal(first[:size].reshape(2,ns),a['initial_canonical_columns']) or np.count_nonzero(first[size:]):
        raise ValueError('fine trajectory did not start from the same upstream columns and zero difference')
    A,D=last[:size].reshape(2,ns),last[size:].reshape(2,ns,N)
    if not np.array_equal(A,a['prepared/reference_amplitudes']):raise ValueError('final reference differs')
    grid,target=a['computational_z'],a['target_z']
    Dt=np.moveaxis(trigonometric_polynomial(D,grid,target),-1,0)
    Dzt=np.moveaxis(trigonometric_polynomial(D,grid,target,derivative=True),-1,0)
    X=A[None]+Dt;phase=np.exp(-1j*a['source/energies']*target[:,None])[:,None,:]
    if not np.array_equal(phase*X,a['prepared/columns']) or not np.array_equal(phase*(Dzt-1j*a['source/energies']*X),a['prepared/axial_columns']):
        raise ValueError('fine incoming fields do not reconstruct from saved trajectory')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    result=run() if p.parse_args().run else check()
    print(json.dumps({k:result[k] for k in ('status','complete_capture','runtime','grid','options')},indent=2))
