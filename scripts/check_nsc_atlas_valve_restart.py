"""Reproduce the bounded same-state local restart to the R≈0.15 cut."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_spherical import (
    SphericalValve,ConformalRestart,align_conformal_radius,conformal_metric_diagnostics,
)

ROOT=Path(__file__).resolve().parents[1]
BASE='226d30e7'
SOURCES=('src/recursive_horizons/nsc_atlas_valve_spherical.py',
         'src/recursive_horizons/nsc_atlas_valve_neck.py',
         'src/recursive_horizons/nsc_atlas_parent_circulation.py',
         'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
         'scripts/check_nsc_atlas_valve_restart.py')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def historical_profile(points=8001):
    """Authenticate the old cache against old Git bytes, not the new module."""
    archived=subprocess.check_output(['git','show',f'{BASE}:results/summaries/nsc-atlas-valve-conformal-v1.json'],cwd=ROOT)
    summary=json.loads(archived);pins=summary['provenance']['source_sha256_after']
    for name,digest in pins.items():
        old=subprocess.check_output(['git','show',f'{BASE}:{name}'],cwd=ROOT)
        if hashlib.sha256(old).hexdigest()!=digest:raise RuntimeError(('historical source pin mismatch',name))
    member=next(v for v in summary['levels'] if v['points']==points)
    cache=ROOT/member['raw_cache']['path']
    if sha(cache)!=member['raw_cache']['sha256']:
        raise RuntimeError('historical profile cache is missing or mismatched; regenerate that historical calculation first')
    with np.load(cache) as raw:
        metadata=json.loads(str(raw['metadata']))
        profile={key:raw[key].copy() for key in ('x','y','yt','yx','yxx')}
        profile.update(T=float(raw['T']),center_proper_time=float(raw['center_proper_time']))
    if metadata['source_sha256']!=pins:raise RuntimeError('historical cache metadata mismatch')
    return profile,dict(git_commit=BASE,archived_summary_sha256=hashlib.sha256(archived).hexdigest(),
                        source_sha256=pins,parent_points=points,raw_cache=member['raw_cache'],
                        rule='Authenticate archived input with historical Git bytes; current analysis sources have separate pins.')


def comparison(values):
    spreads=[abs(values[i+1]-values[i]) for i in (0,1)]
    return dict(values=values,adjacent_spreads=spreads,spread_ratio=spreads[0]/spreads[1] if spreads[1] else None,
                finest_spread_indicator=spreads[1],meaning='Conditional empirical successive-grid difference, not a continuum or propagated-parent-error bound.')


def calculation(cache_directory):
    start=time.perf_counter();before={name:sha(ROOT/name) for name in SOURCES}
    profile,parent=historical_profile()
    with threadpool_limits(limits=1):
        system=SphericalValve(24,epsilon=0.)
        frame=ConformalRestart(profile,H=system.H,kappa=system.kappa,sstar=system.sstar,half_width=.8)
        levels=[];events=[]
        for points in (4001,8001,16001):
            result=frame.evolve(points,endpoint=.3,budget_seconds=60.,relative_geometry_limit=1e-3,
                                relative_source_limit=.02,diagnostic_continue_after_guard=points<16001,radius_floor=.15)
            expected='diagnostic_radius_target' if points<16001 else 'controlled_radius_target'
            if result['status']!=expected:raise RuntimeError((points,result['status'],result['message']))
            event=align_conformal_radius(result,.15);events.append(event)
            directory=cache_directory.resolve();directory.mkdir(parents=True,exist_ok=True)
            path=directory/f'points-{points}.npz'
            np.savez_compressed(path,center_history=result['center_history'],center_columns=np.asarray(result['center_columns']),
                                 **result['end_profile'],metadata=json.dumps(dict(current_analysis_source_sha256=before,
                                                                                historical_parent=parent['raw_cache'],points=points)))
            levels.append(dict(points=points,status=result['status'],final_cut_admitted=result['final_cut_admitted'],
                               first_validity_failure=result['first_validity_failure'],validity_limits=result['validity_limits'],
                               initial=result['rows'][0],last_actual_cut=result['rows'][-1],
                               events=result['events'],radius_0p15=event,
                               raw_cache=dict(path=str(path.relative_to(ROOT)),sha256=sha(path)),
                               elapsed_seconds=result['elapsed_seconds']))
            print(points,result['status'],'R',result['rows'][-1]['center_radius'],
                  'metric/source',result['rows'][-1]['metric_Einstein_residual_over_source_RMS'],flush=True)
        keys=('center_proper_time','center_Rproper','center_Rproper_acc','center_S','center_Pi_S','center_theta_plus_over_H')
        comparisons={key:comparison([v['values'][key] for v in events]) for key in keys}
        # Independent direct-curvature control, distinct from source repair.
        H=system.H;T=.3;x=np.linspace(-.7,.7,1001);zero=np.zeros_like(x)
        scale=1/math.cos(H*T);R=scale*np.cos(H*x)/H;a=-math.log(math.cos(H*T));aT=H*math.tan(H*T)
        y=np.array([R,np.full_like(x,a),zero,zero]);yt=np.array([aT*R,np.full_like(x,aT),zero,zero])
        dx=np.array([-scale*np.sin(H*x),zero,zero,zero]);dxx=np.array([-H*H*R,zero,zero,zero])
        dtt=np.array([(H*H*scale*scale+aT*aT)*R,np.full_like(x,H*H*scale*scale),zero,zero])
        control=conformal_metric_diagnostics(T,0.,y,yt,dx,dxx,x[1]-x[0],H,1.,system.sstar,ytt=dtt)
        control={key:control[key] for key in ('center_metric_Ricci_over_H2','center_metric_Kretschmann_over_H4',
                                            'center_metric_Einstein_residual_over_source','metric_Einstein_residual_over_source_RMS')}
    after={name:sha(ROOT/name) for name in SOURCES}
    if before!=after:raise RuntimeError('analysis sources changed during the bounded calculation')
    return dict(schema='nsc-atlas-valve-restart-v1',historical_input=parent,
                current_analysis=dict(source_sha256_before=before,source_sha256_after=after,source_unchanged=True,
                                      git_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                                      role='New working-tree analysis pins; not authenticated as the historical226d30e7 implementation.'),
                restart=dict(coordinate_scale=frame.scale,half_width=frame.half_width,
                             initial_proper_time=frame.initial_proper_time,
                             rule="T'=c(T-T0),X'=cX,a'=a-log(c); coordinate derivatives scale by1/c or1/c²; fields and physical momenta retained.",
                             representation='C2 cubic interpolation of the actual Cauchy values/velocities; no incoming state or constraint-source replacement.',
                             stored_jet_representation_gaps=frame.representation_gaps),
                levels=levels,radius_0p15_comparison=comparisons,
                direct_de_sitter_control=dict(result=control,exact_Ricci_over_H2=12.,exact_Kretschmann_over_H4=24.),
                invariant_relations=dict(null_product='theta_plus*theta_minus=-4 grad(R)^2/R²',
                                         Hessian='Hess_ss(R)=R_ll+K_rad R_proper',
                                         canonical_NEC='rho+p_r=sum(Pi²+phi_l²)>=0',
                                         turning_point='At R_proper=0 and strict R_ll>0, R_proper_proper=-R_ll-kappa*R*(rho+p_r)/2<0.'),
                numerical=dict(endpoint_radius=.15,new_points=[4001,8001,16001],BLAS_threads=1,
                               temporal_metric_jet='Actual three-time-level differences at evolved cuts. Initial Cauchy temporal jets are PDE-initialized and explicitly not independent.',
                               elapsed_seconds=time.perf_counter()-start,full_arrays_committed=False,
                               scope='Existing characteristic engine and equations; new bounded domain only.'),
                verdict=dict(finest_cut_admitted=levels[-1]['final_cut_admitted'],
                             coarse_radius_events='Diagnostic refinement only: the4k/8k cuts fail the declared source residual guard and are not admitted physical cuts.',
                             finite_errors='Successive new-grid differences are conditional on the trusted numerical parent profile. Its prior finite uncertainty is retained; full propagation of that input uncertainty is not certified.',
                             bounce_observed=any(e['kind']=='center_proper_time_bounce' for v in levels for e in v['events']),
                             neck_type_change_observed=any(e['kind']=='spatial_neck_type_change_on_conformal_slices' for v in levels for e in v['events']),
                             regular_child_proved=False,physical_singularity_proved=False,
                             outcome='The same branch continues regular angular contraction to controlledR≈.15; actual metric curvature grows strongly but remains finite at the stop.',
                             receiver_scope='The local causal boundary is not the whole receiving-band marginal boundary. Zero center flux still permits symmetric feeding; no new late whole-band influx is inferred.'),
                reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_restart.py --output results/summaries/nsc-atlas-valve-restart-v1.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-restart-v1.json')
    parser.add_argument('--cache-directory',type=Path,default=ROOT/'.build/valve-restart-v1')
    args=parser.parse_args();result=calculation(args.cache_directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),radius_event=result['radius_0p15_comparison'],
                          finest_actual_cut=result['levels'][-1]['last_actual_cut']),indent=2))


if __name__=='__main__':main()
