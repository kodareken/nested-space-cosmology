"""Run a bounded nonlinear inherited-parent GH evolution; save a compact record."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from recursive_horizons.nsc_atlas_parent_gh import evolve_parent

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--degree',type=int,default=6);p.add_argument('--epsilon',type=float,default=.025)
    p.add_argument('--proper-endpoint',type=float);p.add_argument('--threads',type=int,default=8)
    p.add_argument('--budget-seconds',type=float,default=600);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--gauge',choices=('harmonic','vacuum_reference'),default='harmonic')
    p.add_argument('--backend',choices=('hopf','dense'),default='hopf')
    p.add_argument('--state-cache',type=Path)
    a=p.parse_args();out=(ROOT/a.output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh repository output required')
    names=['src/recursive_horizons/nsc_atlas_parent_gh.py','src/recursive_horizons/nsc_atlas_s3_frame.py',
        'src/recursive_horizons/nsc_atlas_gh_pointwise.py','src/recursive_horizons/nsc_atlas_gh_constraints.py',
        'src/recursive_horizons/nsc_atlas_finite_parent_data.py','src/recursive_horizons/nsc_atlas_s3_hopf.py',
        'scripts/evolve_nsc_atlas_parent_gh.py']
    before={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    saved=[]
    if a.state_cache:
        cache=(ROOT/a.state_cache).resolve()
        if not cache.is_relative_to(ROOT/'.build') or cache.exists():raise ValueError('fresh state cache under .build required')
    def sink(t,y):saved.append((t,y.copy()))
    r=evolve_parent(degree=a.degree,epsilon=a.epsilon,proper_endpoint=a.proper_endpoint,
        threads=a.threads,budget_seconds=a.budget_seconds,progress=True,gauge=a.gauge,
        state_sink=sink if a.state_cache else None,backend=a.backend)
    after={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    r['source_sha256']=before;r['producer_source_unchanged_during_run']=before==after
    if a.state_cache:
        metadata=dict(degree=a.degree,epsilon=a.epsilon,gauge=a.gauge,backend=a.backend,source_sha256=r['source_sha256'])
        np.savez_compressed(cache,times=np.array([v[0] for v in saved]),states=np.array([v[1] for v in saved]),metadata=json.dumps(metadata))
    out.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=r['status'],stop_reason=r['stop_reason'],degree=r['degree'],epsilon=r['epsilon'],final=r['snapshots'][-1],
        steps=r['steps'],rhs_calls=r['rhs_calls'],elapsed_seconds=r['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
