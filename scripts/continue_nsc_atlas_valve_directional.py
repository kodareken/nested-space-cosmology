"""Carry admitted directional Cauchy data into a shrinking conformal patch.

Preparation changes coordinates only. Evolution uses the package's existing
characteristic diamonds, with reflection explicitly disabled. No incoming
boundary state is supplied at the contracting edges of the causal domain.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import PchipInterpolator
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve,_conformal_evolve

ROOT=Path(__file__).resolve().parents[1]
PHYSICAL='src/recursive_horizons/nsc_atlas_valve_spherical.py'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def physical_symbols(source):
    selected={'canonical_potential','gegenbauer_jets','SphericalValve'}
    nodes=[node for node in ast.parse(source).body if getattr(node,'name',None) in selected]
    if len(nodes)!=len(selected):raise ValueError('spherical decoder symbols missing')
    return [ast.dump(node,include_attributes=False) for node in nodes]


def authentic_run(path):
    path=(ROOT/path).resolve();raw=path.read_bytes();record=json.loads(raw)
    if (not path.is_relative_to(ROOT/'.build') or record['status'] in ('running','initializing') or
        not record['source_bound_to_commit'] or not record['source_unchanged']):
        raise ValueError('terminal source-bound directional run required')
    for name,digest in record['source_sha256'].items():
        old=subprocess.check_output(['git','show',record['producer_commit']+':'+name],cwd=ROOT)
        if hashlib.sha256(old).hexdigest()!=digest:raise ValueError('archived source mismatch: '+name)
        if name.startswith('src/') and sha(ROOT/name)!=digest:
            if name!=PHYSICAL or physical_symbols(old)!=physical_symbols((ROOT/name).read_bytes()):
                raise ValueError('actual physical decoder changed: '+name)
    cache=(ROOT/record['state_cache']).resolve()
    if not cache.is_relative_to(ROOT/'.build') or sha(cache)!=record['state_cache_sha256']:
        raise ValueError('directional Cauchy cache mismatch')
    with np.load(cache,allow_pickle=False) as data:times=data['times'].copy();states=data['states'].copy()
    if states.shape!=(len(times),8,record['modes']) or len(times)!=len(record['snapshots']):
        raise ValueError('aligned Cauchy state history required')
    index=max(i for i,row in enumerate(record['snapshots']) if row['sampled_admitted'])
    row=record['snapshots'][index]
    if abs(times[index]-row['proper_time'])>2e-12:raise ValueError('Cauchy time mismatch')
    return record,row,states[index],dict(path=path.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(raw).hexdigest(),
        state_cache=record['state_cache'],state_cache_sha256=record['state_cache_sha256'],snapshot_index=index,
        proper_time=float(times[index]),scope='Last sampled-admitted state only; no uncached tail is used.')


def source_pins():
    paths={Path(__file__).resolve()};package=(ROOT/'src/recursive_horizons').resolve()
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name:
            path=Path(name).resolve()
            if path.suffix=='.py' and path.is_relative_to(package):paths.add(path)
    return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}


def prepare(path,folder,points,u_max=.9):
    started=time.perf_counter();record,row,state,origin=authentic_run(path)
    system=SphericalValve(record['modes'],reflection_symmetry=False)
    center=row['minimum']['u'] if row['minimum'] else row.get('minimum_anchor_u',0.)
    u=np.linspace(-u_max,u_max,8193);pieces=system.pieces(state,u)
    distance=cumulative_simpson(pieces['q']/np.sqrt(1-u*u),x=u,initial=0.)
    at_center=float(PchipInterpolator(u,distance)(center))
    half_width=.98*min(at_center,float(distance[-1])-at_center)
    signed=at_center-distance
    inverse=PchipInterpolator(signed[::-1],u[::-1])
    X=np.linspace(-half_width,half_width,points);coords=inverse(X);coords[points//2]=center
    z=system.pieces(state,coords);q=z['q'];R=z['radius'];v=1-coords*coords
    y=np.array([R,np.zeros(points),z['y'][4],z['y'][5]])
    yt=np.array([-z['kb']*R,-z['ka'],z['y'][6],z['y'][7]])
    yx=np.array([z['rr']/q,np.zeros(points),-np.sqrt(v)/q*z['yu'][4],-np.sqrt(v)/q*z['yu'][5]])
    yxx=np.array([-z['Rr']*R/2,np.zeros(points),z['lapS']-2*yx[0]*yx[2]/R,z['lapC']-2*yx[0]*yx[3]/R])
    if not np.all(np.isfinite([y,yt,yx,yxx])) or min(R)<=0:raise ValueError('finite positive-radius Cauchy patch required')
    pins=source_pins();raw=folder/'initial.npz'
    np.savez_compressed(raw,x=X,u=coords,y=y,yt=yt,yx=yx,yxx=yxx)
    result=dict(schema='nsc-atlas-directional-conformal-initial-v1',origin=origin,
        source_sha256=pins,git_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        cache=raw.relative_to(ROOT).as_posix(),cache_sha256=sha(raw),points=points,
        center_u=float(center),half_width=half_width,step=float(X[1]-X[0]),
        H=system.H,kappa=system.kappa,sstar=system.sstar,
        physical_mapping='X=proper radial distance increasingchi on the admitted slice; a=0,a_T=-K_rad,R_T=-K_tan*R,S_T=Pi_S,chi_T=Pi_chi.',
        coordinates_only=True,reflection_symmetry=False,incoming_boundary_state_supplied=False,
        minimum_radius=float(min(R)),u_range=[float(min(coords)),float(max(coords))],
        elapsed_seconds=time.perf_counter()-started)
    return result


def evolve(initial_path,folder,endpoint,budget):
    path=(ROOT/initial_path).resolve();initial=json.loads(path.read_bytes());cache=ROOT/initial['cache']
    if initial['schema']!='nsc-atlas-directional-conformal-initial-v1' or sha(cache)!=initial['cache_sha256']:
        raise ValueError('SHA-bound conformal Cauchy jets required')
    for name,digest in initial['source_sha256'].items():
        stored=subprocess.check_output(['git','show',initial['git_parent_commit']+':'+name],cwd=ROOT)
        if hashlib.sha256(stored).hexdigest()!=digest:raise ValueError('Cauchy producer source identity failed: '+name)
    with np.load(cache,allow_pickle=False) as raw:
        grid={k:raw[k].copy() for k in ('y','yt','yx','yxx')}
    grid.update(step=initial['step'],half_width=initial['half_width'])
    pins=source_pins()
    result=_conformal_evolve(grid,H=initial['H'],kappa=initial['kappa'],sstar=initial['sstar'],
        proper_time=initial['origin']['proper_time'],endpoint=endpoint,sample_every=.05,budget_seconds=budget,
        validity_mode='relative_metric_source',relative_geometry_limit=1e-3,relative_source_limit=.02,
        radius_floor=.15,reflection_symmetry=False)
    raw=folder/'history.npz';profile=result.pop('end_profile');history=result.pop('center_history')
    payload=dict(center_history=history)
    if profile is not None:
        payload.update({'profile_'+k:v for k,v in profile.items() if isinstance(v,np.ndarray)})
        result['end_profile_metadata']={k:v for k,v in profile.items() if not isinstance(v,np.ndarray)}
    np.savez_compressed(raw,**payload)
    result.update(schema='nsc-atlas-directional-conformal-evolution-v1',
        initial=path.relative_to(ROOT).as_posix(),initial_sha256=sha(path),initial_origin=initial['origin'],
        source_sha256=pins,source_unchanged=pins=={name:sha(ROOT/name) for name in pins},
        producer_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        raw_cache=raw.relative_to(ROOT).as_posix(),raw_cache_sha256=sha(raw),
        center_clock='Proper time along the new constant-X coordinate worldline; not silently identified with the earlier synchronous geodesic.',
        regular_child_proved=False)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare-run',type=Path);mode.add_argument('--initial',type=Path)
    parser.add_argument('--points',type=int,default=2049);parser.add_argument('--endpoint',type=float,default=.4)
    parser.add_argument('--budget',type=float,default=60.);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();folder=(ROOT/args.output).resolve()
    if not folder.is_relative_to(ROOT/'.build') or folder.exists():raise ValueError('fresh output directory under .build required')
    if args.points<33 or args.points%2!=1:raise ValueError('odd point count>=33 required')
    folder.mkdir(parents=True)
    with threadpool_limits(limits=2):
        result=prepare(args.prepare_run,folder,args.points) if args.prepare_run else evolve(args.initial,folder,args.endpoint,args.budget)
    output=folder/'record.json';output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(output),schema=result['schema'],status=result.get('status'),
        final_row=result.get('rows',[{}])[-1],half_width=result.get('half_width')),indent=2))


if __name__=='__main__':main()
