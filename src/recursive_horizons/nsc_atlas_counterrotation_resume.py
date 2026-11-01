"""Authenticated same-grid terminal checkpoints; no initialization or evolution."""
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np

PRODUCER='scripts/evolve_nsc_atlas_valve_counterrotation.py'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _member(root,name):
    path=(root/name).resolve()
    if not path.is_relative_to(root):raise ValueError('resume member leaves repository')
    return path


def _terminal_allowed(record):
    status=record.get('status');reason=record.get('stop_reason')
    if status=='endpoint_reached':return
    if status!='incomplete' or reason not in (
        'declared formation run budget exhausted','nonlinear GH compute budget exceeded'):
        raise ValueError('resume requires terminal budget stop or completed earlier endpoint; admission failures are not resumable')


@dataclass
class ResumeCheckpoint:
    record: dict
    metadata: dict
    state: np.ndarray
    coordinate_t: float
    proper: float
    baseline_excitation_energy: float
    cumulative_work: float
    provenance: dict


def load_resume_checkpoint(record_path,*,repository_root,degree):
    root=Path(repository_root).resolve();path=Path(record_path).resolve()
    if not path.is_relative_to(root/'.build'):raise ValueError('terminal resume record under .build required')
    raw=path.read_bytes();record=json.loads(raw);_terminal_allowed(record)
    if (record.get('schema')!='nsc-atlas-parent-gh-v3' or record.get('degree')!=degree or
        record.get('source_bound_to_commit') is not True or
        record.get('producer_source_unchanged_during_run') is not True or
        record.get('initial_record_unchanged') is not True or record.get('passive_work_valid') is not True):
        raise ValueError('same-grid terminal source-bound record with valid cumulative work required')
    commit=record.get('producer_commit','');pins=record.get('source_sha256',{})
    required={PRODUCER,'src/recursive_horizons/nsc_atlas_parent_gh.py',
              'src/recursive_horizons/nsc_atlas_gh_pointwise.py',
              'src/recursive_horizons/nsc_atlas_gh_history.py',
              'src/recursive_horizons/nsc_atlas_s3_hopf.py'}
    if not re.fullmatch('[0-9a-f]{40}',commit) or not required.issubset(pins):
        raise ValueError('parent producer and package decoder/RHS pins required')
    for name,digest in pins.items():
        if not re.fullmatch('[0-9a-f]{64}',str(digest)):raise ValueError('invalid parent source hash')
        stored=subprocess.run(['git','show',commit+':'+name],cwd=root,capture_output=True)
        if stored.returncode or hashlib.sha256(stored.stdout).hexdigest()!=digest:
            raise ValueError('parent Git source identity failed: '+name)
        # Only the historical runner may have acquired this resume feature.
        # Every historical package decoder, action and RHS remains identical.
        if name!=PRODUCER and sha(_member(root,name))!=digest:
            raise ValueError('current package/source differs from parent: '+name)
    if sha(_member(root,record['initial_record']))!=record['initial_record_sha256']:
        raise ValueError('original initial preparation record changed')
    cache=_member(root,record['state_cache'])
    if sha(cache)!=record['state_cache_sha256']:raise ValueError('parent sparse cache hash mismatch')
    with np.load(cache,allow_pickle=False) as a:
        times=a['times'].copy();states=a['states'].copy();meta=json.loads(str(a['metadata']))
    snapshots=record.get('snapshots',[])
    if (times.ndim!=1 or states.ndim!=2 or len(times)!=len(states) or len(times)!=len(snapshots) or
        not len(times) or np.any(~np.isfinite(times)) or np.any(~np.isfinite(states)) or np.any(np.diff(times)<=0)):
        raise ValueError('invalid aligned sparse checkpoint arrays')
    if meta.get('source_sha256')!=pins or meta.get('producer_commit')!=commit:
        raise ValueError('sparse decoder identities differ from parent record')
    for key,expected in (('degree',degree),('epsilon',record['epsilon']),('gauge',record['prescribed_GH_gauge']),
        ('backend',record['spatial_transform']),('kick_amplitude',record['model_inputs']['kick_amplitude']),('kick_source',2),
        ('preparation_id','one-mode-counterrotation-v1')):
        if meta.get(key)!=expected:raise ValueError('same preparation/decoder required: '+key)
    for t,state,snapshot in zip(times,states,snapshots):
        if (not math.isclose(float(t),snapshot['coordinate_t'],rel_tol=0,abs_tol=2e-12) or
            not math.isclose(float(state[-2]),snapshot['observer_proper_duration'],rel_tol=0,abs_tol=2e-10)):
            raise ValueError('sparse state does not match diagnostic event')
    admitted=[i for i,s in enumerate(snapshots) if s.get('sampled_admitted') is True]
    if not admitted:raise ValueError('no admitted sparse checkpoint')
    index=admitted[-1];snapshot=snapshots[index];t=float(times[index]);state=states[index].copy()
    if index!=len(snapshots)-1 or not math.isclose(t,record['last_admitted_coordinate_t'],rel_tol=0,abs_tol=2e-12):
        raise ValueError('terminal record has an unadmitted later diagnostic; do not resume same grid')
    baseline=float(record.get('excitation_work_origin_energy',snapshots[0]['normal_energy_minus_fixed_vacuum']))
    work=float(snapshot['integrated_excitation_work'])
    if not math.isfinite(baseline+work):raise ValueError('finite original energy baseline/cumulative work required')
    prefix=[]
    for item in record.get('physical_history_prefix',[])+record.get('dense_history_chunks',[]):
        chunk=_member(root,item['path'])
        if sha(chunk)!=item['sha256']:raise ValueError('parent dense history hash failed')
        with np.load(chunk,allow_pickle=False) as a:
            bounds=a['bounds'];chunkmeta=json.loads(str(a['metadata']))
        if bounds.ndim!=2 or bounds.shape[1]!=2 or not len(bounds):raise ValueError('invalid parent history bounds')
        declared=item.get('usable_coordinate_bounds',item['coordinate_bounds'])
        lo=max(float(bounds[0,0]),float(declared[0]));hi=min(float(bounds[-1,1]),float(declared[1]),t)
        if hi>lo:prefix.append(dict(path=item['path'],sha256=item['sha256'],
            coordinate_bounds=[float(bounds[0,0]),float(bounds[-1,1])],usable_coordinate_bounds=[lo,hi],
            producer_commit=chunkmeta['producer_commit']))
    if not prefix or prefix[0]['usable_coordinate_bounds'][0]!=0. or any(
        b['usable_coordinate_bounds'][0]!=a['usable_coordinate_bounds'][1] for a,b in zip(prefix,prefix[1:])) or prefix[-1]['usable_coordinate_bounds'][1]!=t:
        raise ValueError('inherited physical history must continuously cover only the original start through handoff')
    provenance=dict(parent_run=path.relative_to(root).as_posix(),parent_run_sha256=hashlib.sha256(raw).hexdigest(),
        parent_producer_commit=commit,parent_status=record['status'],parent_stop_reason=record.get('stop_reason'),
        sparse_snapshot_index=index,handoff_coordinate_t=t,handoff_proper=float(state[-2]),
        state_array_sha256=hashlib.sha256(state.tobytes()).hexdigest(),physical_history_prefix=prefix,
        preserved_unpromoted_parent_tail=dict(coordinate_bounds=[t,record['actual_last_accepted_coordinate_t']],
            accepted_endpoint_cache=record.get('accepted_endpoint_cache'),
            scope='Original tail files preserved; excluded after handoff, never relabelled as admitted.'),
        baseline_excitation_energy=baseline,cumulative_work_at_handoff=work)
    return ResumeCheckpoint(record,meta,state,t,float(state[-2]),baseline,work,provenance)


def validate_resumed_system(checkpoint,system):
    if checkpoint.state.shape!=(system.dof+10,):raise ValueError('resume full packed state/observer shape mismatch')
    model=checkpoint.record['model_inputs']
    for key,actual in (('H_reference',system.r),('vacuum',system.parent.vacuum),
        ('scalar_mass_squared',system.parent.mass_squared),('volume_kinetic',system.K),('sstar',system.sstar)):
        if not math.isclose(model[key],actual,rel_tol=2e-13,abs_tol=2e-14):raise ValueError('resume physical model changed: '+key)


def remaining_proper_targets(endpoint,handoff_proper):
    if endpoint<=handoff_proper:raise ValueError('resume target must be later than the admitted checkpoint')
    targets=sorted(set([.02]+[i/10 for i in range(1,math.ceil(10*endpoint)+1)]+[endpoint]))
    return [v for v in targets if handoff_proper+1e-10<v<=endpoint]
