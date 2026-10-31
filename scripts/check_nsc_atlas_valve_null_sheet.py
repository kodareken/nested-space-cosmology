"""Trace an inherited closed null sheet on the actual saved parent geometry.

No new GH evolution or Einstein-source curvature replacement is performed.
Both expansions refer to each same common-time cut. Caustics and finite
data boundaries are recorded without inventing a receiving child.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_null_parent import GHNullHistorySampler
from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_null_cuts import initial_null_bundle,evaluate_parametric_cut
from recursive_horizons.nsc_atlas_null_transport import transport_null_bundle
from recursive_horizons.nsc_atlas_trapped_surfaces import StarShapedSurface,round_surface,_unit_frame

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def coordinate_at_proper(sampler,proper):
    history=sampler.history;dof=sampler.system.dof
    if proper==0:return 0.
    endpoint=sampler.data_end
    value=float(history.evaluate(endpoint)[dof+8])
    if abs(value-proper)<2e-10:return endpoint
    return brentq(lambda t:float(history.evaluate(t)[dof+8])-proper,0.,endpoint,xtol=1e-10)


def protected_plane_probes(sampler):
    """Four exact nodal generators on the reflection-fixed initial equator."""
    X=np.array([[0.,0.,0.,1.],[0.,0.,0.,-1.],[0.,1.,0.,0.],[0.,-1.,0.,0.]])
    d=sampler.slice(0.,X);h=d['h'];hi=np.linalg.inv(h);frame=_unit_frame(X)
    # Outward from center X3=+1 at rho=pi/2 is minus the X3 gradient.
    cov=-sampler.r*frame[:,:,2]
    normal=np.einsum('nij,nj->ni',hi,cov)
    normal/=np.sqrt(np.einsum('ni,ni->n',cov,normal))[:,None]
    screens=np.empty((4,2,3));P=np.empty((4,2,2))
    for i in range(4):
        L=np.linalg.cholesky(h[i]);unit=L.T@normal[i]
        axis=np.eye(3)[int(np.argmin(abs(unit)))];e=axis-unit*(axis@unit);e/=np.linalg.norm(e)
        screens[i]=np.linalg.solve(L.T,np.stack((e,np.cross(unit,e)),axis=1)).T
        F=float(np.trace(h[i])/3)
        if np.max(abs(h[i]-F*np.eye(3)))>1e-10*F:
            raise ValueError('protected initial-plane probes require the declared conformal initial geometry')
        derivative=np.trace(d['dh'][i],axis1=1,axis2=2)/3
        b=float(normal[i]@derivative)/(2*F)
        P[i]=b*np.eye(2)-screens[i]@d['K'][i]@screens[i].T
    n=np.c_[1/d['lapse'],-d['shift']/d['lapse'][:,None]]
    return dict(X=X,k=n+np.c_[np.zeros(4),normal],screen=screens,
        D=np.broadcast_to(np.eye(2),(4,2,2)).copy(),P=P)


def initial_sheet(sampler,seed,order):
    if seed=='reflection':
        surface=round_surface([0,0,1,0],math.pi/2);t0=0.
    else:
        record=json.loads((ROOT/'results/summaries/nsc-atlas-valve-formation-early-v1.json').read_bytes())
        row=next(r for r in record['surface_searches'] if r['center_name']=='scalar_anchor')
        c=row['selected_candidate']
        surface=StarShapedSurface(row['center'],c['shape_coefficients'],c['shape_degree'],c['radius_bounds'])
        t0=coordinate_at_proper(sampler,.225)
    bundle=initial_null_bundle(surface,lambda X:sampler.slice(t0,X),frame_scale=sampler.r,order=order)
    initial=bundle['outgoing'];angular_rays=len(initial['X']);probes=0
    if seed=='reflection':
        extra=protected_plane_probes(sampler);probes=len(extra['X'])
        initial={name:np.concatenate([initial[name],extra[name]],axis=0) for name in ('X','k','screen','D','P')}
    return t0,bundle,initial,angular_rays,probes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path,required=True)
    parser.add_argument('--seed',choices=('reflection','marginal_anchor'),default='reflection')
    parser.add_argument('--angular-order',type=int,default=10)
    parser.add_argument('--cut-degrees',type=int,nargs='+',default=[6,8])
    parser.add_argument('--budget-seconds',type=float,default=600.)
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--max-step',type=float,default=.03)
    parser.add_argument('--rtol',type=float,default=3e-8)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--ray-cache',type=Path)
    parser.add_argument('--ray-dense-cache',type=Path)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    if args.angular_order<3 or max(args.cut_degrees)>=args.angular_order:
        raise ValueError('angular samples must resolve the cut harmonics')
    if args.seed=='marginal_anchor' and args.angular_order<9:
        raise ValueError('marginal seed shape requires at least order9')
    if not 0<args.threads<=8 or not math.isfinite(args.budget_seconds) or args.budget_seconds<=0:
        raise ValueError('bounded positive computation controls required')
    names=['scripts/check_nsc_atlas_valve_null_sheet.py',
        'src/recursive_horizons/nsc_atlas_null_geometry.py','src/recursive_horizons/nsc_atlas_null_parent.py',
        'src/recursive_horizons/nsc_atlas_null_cuts.py','src/recursive_horizons/nsc_atlas_null_transport.py',
        'src/recursive_horizons/nsc_atlas_trapped_parent.py','src/recursive_horizons/nsc_atlas_trapped_surfaces.py']
    before={name:sha(ROOT/name) for name in names};commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    bound=all(hashlib.sha256(subprocess.check_output(['git','show',commit+':'+name])).hexdigest()==expected
        for name,expected in before.items())
    started=time.perf_counter()
    with threadpool_limits(limits=args.threads):
        sampler=GHNullHistorySampler.from_files(ROOT/args.record)
        order=(args.angular_order,2*args.angular_order+1)
        t0,bundle,initial,angular_rays,probes=initial_sheet(sampler,args.seed,order)
        proper_points=[p for p in (.075,.15,.225,.3,.6,.9,1.2) if coordinate_at_proper(sampler,p)>t0+1e-10]
        checkpoints=[coordinate_at_proper(sampler,p) for p in proper_points]
        print('Tracing',len(initial['X']),'rays from t=',t0,'to',sampler.data_end,flush=True)
        ray_history=DenseGHHistory() if args.ray_dense_cache else None
        transport=transport_null_bundle(sampler,initial,coordinate_interval=(t0,sampler.data_end),
            frame_scale=sampler.r,checkpoints=checkpoints,rtol=args.rtol,atol=args.rtol/100,
            max_step=args.max_step,budget_seconds=args.budget_seconds,progress=True,data_end=sampler.data_end,
            dense_sink=None if ray_history is None else ray_history.record)
        result=dict(schema='nsc-atlas-valve-null-sheet-v1',seed=args.seed,
            metric_degree=sampler.system.basis.degree,epsilon=sampler.system.epsilon,
            producer_commit=commit,source_sha256=before,source_bound_to_commit=bound,
            parent_history=sampler.provenance,angular_order=list(order),cut_degrees=args.cut_degrees,
            angular_rays=angular_rays,protected_probe_rays=probes,
            transport=transport.summary(),cuts=[],physical_child_admitted=False,
            scope='A finite transported null sheet on stored parent geometry. Sampled cut signs and probe obstructions are not a global no-formation or a child-continuation theorem.')
        for t,state,diagnostic in zip(transport.times,transport.states,transport.diagnostics):
            row=dict(coordinate_t=float(t),cuts=[],normalization='theta_affine=(N*k0)*theta_unit_normal')
            data=sampler(t,state['X']);k=state['k'];energy=data['lapse']*k[:,0]
            source=np.einsum('nab,na,nb->n',data['matter_Ricci'],k,k)
            actual=np.einsum('nab,na,nb->n',data['Ricci'],k,k)
            row.update(metric_Ricci_kk_range=[float(actual.min()),float(actual.max())],
                scalar_source_Ricci_kk_range=[float(source.min()),float(source.max())],
                max_abs_null_Ricci_source_gap=float(np.max(abs(actual-source))),
                max_metric_velocity_consistency=float(np.max(abs(data['metric_velocity_consistency']))))
            if probes:
                row['protected_probes']=dict(points=state['X'][-probes:].tolist(),
                    affine_parameter=state['affine_lambda'][-probes:].tolist(),
                    chi=data['q'][-probes:,1].tolist(),Pi_chi=data['Pi'][-probes:,1].tolist(),
                    Pi_S=data['Pi'][-probes:,0].tolist(),
                    tangential_capacity_over_H=(data['future_tangential_capacity'][-probes:]/sampler.r).tolist(),
                    normal_expansion_over_H=(-np.sum(data['K_eigenvalues'][-probes:],axis=1)/sampler.r).tolist(),
                    unit_null_expansion=[float(x) if np.isfinite(x) else None
                        for x in diagnostic['theta_over_normal_energy'][-probes:]],
                    scope='These generator events obstruct cuts through those events if capacity is negative; no claim about all affine cuts or other sheets.')
            caustic=bool(np.any(diagnostic['caustic_rays_so_far']))
            for degree in args.cut_degrees:
                if float(t)==t0:
                    row['cuts'].append(dict(angular_degree=degree,initial_exact=True,
                        theta_plus_range=[float(bundle['theta_plus'].min()),float(bundle['theta_plus'].max())],
                        theta_minus_range=[float(bundle['theta_minus'].min()),float(bundle['theta_minus'].max())],
                        area=bundle['area'],physical_child_admitted=False))
                    continue
                try:
                    cut=evaluate_parametric_cut(state['X'][:angular_rays],lambda X:sampler.slice(float(t),X),
                        frame_scale=sampler.r,angular_degree=degree,order=order,
                        outgoing_k=state['k'][:angular_rays],caustic_encountered=caustic,
                        ray_area_determinant=np.linalg.det(state['D'][:angular_rays]),expansion_tolerance=1e-6*sampler.r)
                    summary=cut.summary();theta=diagnostic['theta_over_normal_energy'][:angular_rays]
                    good=diagnostic['theta_defined'][:angular_rays]
                    if np.all(good):
                        gap=cut.theta_plus-theta
                        summary['Sachs_geometric_theta_RMS_gap_over_H']=float(np.sqrt(np.mean(gap**2))/sampler.r)
                        summary['Sachs_geometric_theta_max_gap_over_H']=float(np.max(abs(gap))/sampler.r)
                    summary['Sachs_area']=float(bundle['area_weights']@np.linalg.det(state['D'][:angular_rays]))
                    summary['geometric_vs_Sachs_area_relative_gap']=summary['area']/summary['Sachs_area']-1
                    row['cuts'].append(summary)
                except ValueError as error:row['cuts'].append(dict(angular_degree=degree,rejected_reason=str(error),physical_child_admitted=False))
            result['cuts'].append(row)
        if args.ray_cache:
            cache=(ROOT/args.ray_cache).resolve()
            if not cache.is_relative_to(ROOT/'.build') or cache.exists():raise ValueError('fresh local ray cache required')
            np.savez_compressed(cache,times=transport.times,states=transport.packed_states)
            result['ray_cache']=dict(path=cache.relative_to(ROOT).as_posix(),sha256=sha(cache))
        if ray_history is not None:
            cache=(ROOT/args.ray_dense_cache).resolve()
            if not cache.is_relative_to(ROOT/'.build') or cache.exists():raise ValueError('fresh local dense ray cache required')
            ray_history.save(cache,dict(schema='nsc-null-ray-dense-v1',seed=args.seed,
                rays=len(initial['X']),angular_rays=angular_rays,protected_probe_rays=probes,
                angular_order=list(order),source_sha256=before,producer_commit=commit,
                parent_history=sampler.provenance,coordinate_interval=transport.effective_interval,
                initial_area_weights=bundle['area_weights'].tolist(),
                initial_area=bundle['area']))
            result['ray_dense_cache']=dict(path=cache.relative_to(ROOT).as_posix(),sha256=sha(cache),
                intervals=len(ray_history.bounds),schema='nsc-null-ray-dense-v1')
        result.update(elapsed_seconds=time.perf_counter()-started,sampler_calls=sampler.calls,
            source_unchanged_during_run=before=={name:sha(ROOT/name) for name in names},
            sampler_points=sampler.points_evaluated,sampler_elapsed_seconds=sampler.elapsed_seconds)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print('Saved',output,'transport=',transport.status,'in',result['elapsed_seconds'],'seconds',flush=True)


if __name__=='__main__':main()
