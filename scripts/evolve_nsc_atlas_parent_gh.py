"""Run a bounded nonlinear inherited-parent GH evolution; save a compact record."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np

from recursive_horizons.nsc_atlas_parent_gh import evolve_parent
from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--degree',type=int,default=6);p.add_argument('--epsilon',type=float,default=.025)
    p.add_argument('--proper-endpoint',type=float);p.add_argument('--threads',type=int,default=8)
    p.add_argument('--budget-seconds',type=float,default=600);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--gauge',choices=('harmonic','vacuum_reference','lapse_controlled'),default='lapse_controlled')
    p.add_argument('--backend',choices=('hopf','dense'),default='hopf')
    p.add_argument('--state-cache',type=Path)
    p.add_argument('--dense-cache',type=Path)
    p.add_argument('--kick-amplitude',type=float,default=0.)
    p.add_argument('--kick-source',type=int,choices=(1,2),default=1)
    p.add_argument('--quadrature-radial',type=int);p.add_argument('--quadrature-angles',type=int)
    p.add_argument('--rtol',type=float,default=3e-8);p.add_argument('--atol',type=float,default=3e-10)
    p.add_argument('--max-step',type=float,default=.05)
    p.add_argument('--no-work',action='store_true',
                   help='Omit passive work quadrature for the formation-only campaign; physical equations are unchanged.')
    a=p.parse_args();out=(ROOT/a.output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh repository output required')
    names=['src/recursive_horizons/nsc_atlas_parent_gh.py','src/recursive_horizons/nsc_atlas_s3_frame.py',
        'src/recursive_horizons/nsc_atlas_gh_pointwise.py','src/recursive_horizons/nsc_atlas_gh_constraints.py',
        'src/recursive_horizons/nsc_atlas_finite_parent_data.py','src/recursive_horizons/nsc_atlas_s3_hopf.py',
        'src/recursive_horizons/nsc_atlas_gh_observer.py','src/recursive_horizons/nsc_atlas_parent_circulation.py',
        'src/recursive_horizons/nsc_atlas_parent_scale_selector.py','src/recursive_horizons/nsc_atlas_parent_backreaction.py',
        'src/recursive_horizons/nsc_atlas_s3_spectral.py',
        'src/recursive_horizons/nsc_atlas_gh_matter_work.py','src/recursive_horizons/nsc_atlas_passive_dop853.py',
        'src/recursive_horizons/nsc_atlas_gh_history.py',
        'src/recursive_horizons/nsc_atlas_kick_york.py','src/recursive_horizons/nsc_atlas_gravitating_kick.py',
        'scripts/evolve_nsc_atlas_parent_gh.py']
    before={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bound=True
    for name in names:
        stored=subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True)
        bound=bound and stored.returncode==0 and hashlib.sha256(stored.stdout).hexdigest()==before[name]
    saved=[]
    history=None
    if a.dense_cache:
        dense_cache=(ROOT/a.dense_cache).resolve()
        if not dense_cache.is_relative_to(ROOT/'.build') or dense_cache.exists():raise ValueError('fresh dense cache under .build required')
        history=DenseGHHistory()
    if a.state_cache:
        cache=(ROOT/a.state_cache).resolve()
        if not cache.is_relative_to(ROOT/'.build') or cache.exists():raise ValueError('fresh state cache under .build required')
    def sink(t,y):saved.append((t,y.copy()))
    if (a.quadrature_radial is None)!=(a.quadrature_angles is None):raise ValueError('both quadrature orders are required')
    order=None if a.quadrature_radial is None else (a.quadrature_radial,a.quadrature_angles)
    r=evolve_parent(degree=a.degree,epsilon=a.epsilon,proper_endpoint=a.proper_endpoint,
        threads=a.threads,budget_seconds=a.budget_seconds,progress=True,gauge=a.gauge,
        state_sink=sink if a.state_cache else None,backend=a.backend,quadrature_order=order,
        rtol=a.rtol,atol=a.atol,max_step=a.max_step,dense_sink=None if history is None else history.record,
        kick_amplitude=a.kick_amplitude,kick_source=a.kick_source,measure_work=not a.no_work)
    after={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    r['source_sha256']=before;r['producer_source_unchanged_during_run']=before==after
    r['producer_commit']=commit;r['source_bound_to_commit']=bound
    if history is not None:
        metadata=dict(degree=a.degree,epsilon=a.epsilon,gauge=a.gauge,backend=a.backend,
            kick_amplitude=a.kick_amplitude,kick_source=a.kick_source,
            quadrature_order=order,temporal_tolerances=r['temporal_tolerances'],
            source_sha256=before,producer_commit=commit,
            physical_coordinate_endpoint=r['snapshots'][-1]['coordinate_t'],
            proper_endpoint=r['snapshots'][-1]['observer_proper_duration'])
        history.save(dense_cache,metadata)
        r['dense_history']=dict(schema='nsc-atlas-gh-dense-v1',intervals=len(history.bounds),
            state_size=len(history.old[0]),sha256=hashlib.sha256(dense_cache.read_bytes()).hexdigest(),
            scope='Local regenerable dense polynomials; no sparse-snapshot interpolation or changed adaptive state.')
    if a.state_cache:
        metadata=dict(degree=a.degree,epsilon=a.epsilon,gauge=a.gauge,backend=a.backend,
            kick_amplitude=a.kick_amplitude,kick_source=a.kick_source,source_sha256=r['source_sha256'])
        np.savez_compressed(cache,times=np.array([v[0] for v in saved]),states=np.array([v[1] for v in saved]),metadata=json.dumps(metadata))
    out.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=r['status'],stop_reason=r['stop_reason'],degree=r['degree'],epsilon=r['epsilon'],final=r['snapshots'][-1],
        steps=r['steps'],rhs_calls=r['rhs_calls'],elapsed_seconds=r['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
