"""Bounded closed-surface formation diagnostics on authenticated saved parents.

No GH evolution is performed. Surface searches use the finest supplied
metric; the same resulting surfaces are evaluated on every supplied metric.
Misses exclude neither other surfaces nor later formation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import time

import numpy as np
from numba import set_num_threads
from scipy.fft import set_workers
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_backreaction import anchor_point
from recursive_horizons.nsc_atlas_parent_gh import FRAME_CHAR,PAIRS
from recursive_horizons.nsc_atlas_trapped_parent import (
    GHSurfaceGeometry,_first_geometry,compare_surface_neighborhood,
)
from recursive_horizons.nsc_atlas_trapped_surfaces import (
    StarShapedSurface,evaluate_surface,find_mots,round_surface,
)

ROOT=Path(__file__).resolve().parents[1]
START_RADII=(.7,1.4,2.0)


def stabilization_threshold(parent):
    """Same-action local barrier loss; not a child or escape criterion."""
    from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
    model=VolumeParent()
    minimum=next(row for row in model.vacuum_roots() if row['kind']=='minimum')
    sstar=minimum['s']
    x=brentq(lambda x:225*x**6-560*x+336,.9,1.05,xtol=1e-14)
    sfold=-math.log(x)/2
    critical=-14*x**3+17.5*x**4-3.125*x**9
    chi_squared_critical=((175/24)*(14/15)**.2-7)/200
    z=parent.mass_squared/parent.H**2
    a=math.sqrt(12/(parent.kappa*(z+3)))
    b=math.sqrt(72/(parent.kappa*(z+8)))
    maximum=(a*a+b*b)/(2*b)
    local_loading=parent.mass_squared*maximum**2/2*math.exp(6*sstar)
    return dict(fold_x=x,fold_s=sfold,critical_adiabatic_loading=critical,
        maximum_initial_field_at_epsilon1=maximum,
        maximum_initial_adiabatic_loading_at_epsilon1=local_loading,
        critical_epsilon_at_scalar_maximum=math.sqrt(critical/local_loading),
        epsilon1_loading_over_critical=local_loading/critical,
        instantaneous_chi_squared_critical=chi_squared_critical,
        instantaneous_absolute_chi_critical=math.sqrt(chi_squared_critical),
        instantaneous_critical_epsilon_at_scalar_maximum=math.sqrt(chi_squared_critical)/maximum,
        E0_reference=3/8,E0_relative_below_critical_percent=100*(1-(3/8)/critical),
        derivation='Veff=V+M*exp(-6s); fold225x^6-560x+336=0,x=exp(-2s). M(s)=-14x^3+17.5x^4-3.125x^9.',
        scope='Fixed local oscillator action and fixed external volume, neglecting spatial transport. Conditional adiabatic approximation, not dynamical escape, a trapped surface or received child energy.')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Budget:
    def __init__(self,seconds):self.started=time.perf_counter();self.seconds=seconds
    @property
    def elapsed(self):return time.perf_counter()-self.started
    def check(self):
        if self.elapsed>=self.seconds:raise TimeoutError('formation diagnostic wall budget exhausted')


class LimitedGeometry:
    """One-point-set cache and a budget guard; never changes the saved state."""
    def __init__(self,provider,budget):
        self.provider=provider;self.budget=budget;self.r=provider.r
        self.points=None;self.data=None;self.calls=0
    def __call__(self,points):
        self.budget.check();X=np.asarray(points)
        if self.points is not None and np.array_equal(X,self.points):return self.data
        self.calls+=1;data=self.provider(X);self.budget.check()
        self.points=X.copy();self.data=data
        return data


def _focusing(data,kappa):
    W=np.asarray(data['W']);kinetic=np.sum(data['Pi']**2,axis=1)
    K2=np.sum(data['K_eigenvalues']**2,axis=1)
    if not all(np.all(np.isfinite(a)) for a in (W,kinetic,K2)):
        raise ValueError('derived focusing inputs are not finite')
    if np.any(W<=0):return None,None
    return (kinetic+K2/kappa)/W,kinetic/W


def grid_diagnostic(system,state,t,index,producer_point,budget):
    budget.check();system.update_point(state[system.dof:system.dof+4])
    C=system.unpack(state);V=system.evaluate(C,system.chars)
    first=np.stack([system.evaluate(system.r*(D@C.T).T,system.chars*FRAME_CHAR[i])
        for i,D in enumerate(system.basis.generators)],axis=1)
    count=len(V);g=np.empty((count,4,4));v=np.empty_like(g);dg=np.empty((count,3,4,4))
    for k,(a,b) in enumerate(PAIRS):
        g[:,a,b]=g[:,b,a]=V[:,k];v[:,a,b]=v[:,b,a]=V[:,k+10]
        dg[:,:,a,b]=dg[:,:,b,a]=first[:,:,k]
    geo=_first_geometry(g,v,dg,system.r)
    geo['Pi']=(V[:,22:24]-np.einsum('ni,nia->na',geo['shift'],first[:,:,20:22]))/geo['lapse'][:,None]
    geo['W']=system.potential(V[:,20:22]);F,Fkin=_focusing(geo,system.parent.kappa)
    s=system.sstar+V[:,20]/math.sqrt(system.K)
    mu2=400*np.exp(-12*s)
    oscillator_loading=(geo['Pi'][:,1]**2+mu2*V[:,21]**2)/2*np.exp(6*s)
    threshold=stabilization_threshold(system.parent)
    loading_ratio=oscillator_loading/threshold['critical_adiabatic_loading']
    volume_weights=system.basis.weights*np.sqrt(np.linalg.det(geo['h'][:-1]))
    theta=-np.sum(geo['K_eigenvalues'],axis=1);cap=geo['future_tangential_capacity'];r=system.r
    constraints={key:producer_point[key] for key in ('ADM_H_RMS_over_H_squared',
        'ADM_M_RMS_over_H_squared','GH_C_RMS_over_H','metric_valid') if key in producer_point}
    row=dict(snapshot=index,coordinate_t=float(t),observer_proper_duration=float(state[system.dof+8]),
        grid_points=system.N,W_range=[float(geo['W'][:-1].min()),float(geo['W'][:-1].max())],
        theta_over_H_reference_range=[float(theta[:-1].min()/r),float(theta[:-1].max()/r)],
        future_tangential_capacity_over_H_reference_range=[float(cap[:-1].min()/r),float(cap[:-1].max()/r)],
        proper_volume_times_H_reference_cubed=float(volume_weights.sum()),producer_constraints=constraints,
        F_defined=F is not None,F_range=None,F_observer=None,F_fixed_scalar_anchor=None,
        kinetic_over_W_max=None,proper_volume_fraction_F_gt_one=None,max_F_point=None)
    row.update(volume_mode_range=[float(s[:-1].min()),float(s[:-1].max())],
        observer_internal_radius_over_initial=math.exp(s[-1]-system.sstar),
        observer_internal_volume_over_initial=math.exp(12*(s[-1]-system.sstar)),
        adiabatic_loading_over_critical_range=[float(loading_ratio[:-1].min()),float(loading_ratio[:-1].max())],
        adiabatic_loading_over_critical_at_observer=float(loading_ratio[-1]))
    normal_s=geo['Pi'][:,0]/math.sqrt(system.K)
    external_expansion=np.exp(6*s)*(theta-18*normal_s)
    row.update(observer_normal_volume_roll=float(normal_s[-1]),
        observer_Einstein_three_expansion_over_H=float(theta[-1]/r),
        observer_external_Jordan_three_expansion_over_H=float(external_expansion[-1]/r),
        external_Jordan_three_expansion_over_H_range=[float(external_expansion[:-1].min()/r),float(external_expansion[:-1].max()/r)],
        frame_scope='Declared block-diagonal d12 breathing-mode uplift. Internal-volume growth is not a newly formed regional child.')
    if F is not None:
        row.update(F_range=[float(F[:-1].min()),float(F[:-1].max())],F_observer=float(F[-1]),
            kinetic_over_W_max=float(Fkin[:-1].max()),
            proper_volume_fraction_F_gt_one=float(volume_weights@(F[:-1]>1)/volume_weights.sum()),
            max_F_point=system.basis.points[int(np.argmax(F[:-1]))].tolist())
        point=GHSurfaceGeometry(system,state,coordinate_t=float(t))(anchor_point(system.parent)[None,:])
        pointF,_=_focusing(point,system.parent.kappa)
        row['F_fixed_scalar_anchor']=None if pointF is None else float(pointF[0])
    budget.check();return row


def surface_diagnostic(evaluation,provider):
    data=provider(evaluation.points);inverse=np.linalg.inv(data['h'])
    kinetic=np.sum(data['Pi']**2,axis=1)
    gradient=np.einsum('nij,nia,nja->n',inverse,data['dq'],data['dq'])
    rho=(kinetic+gradient)/2+data['W'];pressure=(kinetic-gradient/3)/2-data['W']
    momentum=-np.einsum('na,nia->ni',data['Pi'],data['dq'])
    flux=np.einsum('ni,ni->n',momentum,evaluation.outward_normal)
    weights=evaluation.area_weights/evaluation.area
    system=provider.provider.system
    s=system.sstar+data['q'][:,0]/math.sqrt(system.K)
    ns=data['Pi'][:,0]/math.sqrt(system.K)
    vs=np.einsum('ni,ni->n',data['dq'][:,:,0],evaluation.outward_normal)/math.sqrt(system.K)
    factor=np.exp(6*s)
    wrapped_plus=factor*evaluation.theta_plus;wrapped_minus=factor*evaluation.theta_minus
    external_plus=factor*(evaluation.theta_plus-12*(ns+vs))
    external_minus=factor*(evaluation.theta_minus-12*(ns-vs))
    return dict(normal_flux_range=[float(flux.min()),float(flux.max())],
        outward_normal_energy_flux_area_average=float(weights@flux),
        rho_area_average=float(weights@rho),isotropic_pressure_area_average=float(weights@pressure),
        W_area_average=float(weights@data['W']),
        wrapped_16_theta_plus_range=[float(wrapped_plus.min()),float(wrapped_plus.max())],
        wrapped_16_theta_minus_range=[float(wrapped_minus.min()),float(wrapped_minus.max())],
        external_Jordan_only_theta_plus_range=[float(external_plus.min()),float(external_plus.max())],
        external_Jordan_only_theta_minus_range=[float(external_minus.min()),float(external_minus.max())],
        wrapped_null_signs_equal_Einstein=True,
        scope='Instantaneous field stress on a sampled surface; flux positive outward. No child energy, thermal state or vacuum fraction is identified.')


def classification(evaluation):
    if evaluation.plus_max_abs>evaluation.expansion_tolerance:return 'not_converged'
    if evaluation.future_mots_sampled:return 'future_MOTS_sampled'
    return 'marginal_root_not_future_MOTS'


def embed_degree(surface,degree):
    if degree<surface.degree:raise ValueError('shape warm start cannot discard coefficients')
    coefficients=np.zeros((degree+1)**2);coefficients[:len(surface.coefficients)]=surface.coefficients
    return StarShapedSurface(surface.center,coefficients,degree,surface.radius_bounds)


def candidate_summary(candidate):
    surface,evaluation,origin=candidate
    return dict(evaluation.summary(),origin=origin,classification=classification(evaluation),
        shape_degree=surface.degree,shape_coefficients=surface.coefficients.tolist(),
        radius_bounds=list(surface.radius_bounds))


def compare_candidate(surface,providers,degrees,order,tolerance,budget):
    budget.check();evaluations=[evaluate_surface(surface,p,frame_scale=p.r,quadrature_order=order,
        expansion_tolerance=tolerance) for p in providers]
    fine=evaluations[-1];weights=fine.area_weights/fine.area
    tensors=compare_surface_neighborhood(providers,fine.points,weights=fine.area_weights)
    geometry=[p(fine.points) for p in providers];gaps=[]
    for index,(old,new) in enumerate(zip(evaluations,evaluations[1:])):
        g0,g1=geometry[index:index+2];r=providers[index+1].r;kappa=providers[index+1].provider.system.parent.kappa
        F0,_=_focusing(g0,kappa);F1,_=_focusing(g1,kappa)
        plus=new.theta_plus-old.theta_plus;minus=new.theta_minus-old.theta_minus
        gaps.append(dict(from_metric_degree=degrees[index],to_metric_degree=degrees[index+1],
            theta_plus_RMS_gap_over_H_reference=float(math.sqrt(weights@(plus*plus))/r),
            theta_plus_max_abs_gap_over_H_reference=float(max(abs(plus))/r),
            theta_minus_RMS_gap_over_H_reference=float(math.sqrt(weights@(minus*minus))/r),
            theta_minus_max_abs_gap_over_H_reference=float(max(abs(minus))/r),
            canonical_field_RMS_gap_times_sqrt_kappa=np.sqrt(kappa*(weights@((g1['q']-g0['q'])**2))).tolist(),
            normal_momentum_RMS_gap_times_sqrt_kappa_over_H=(np.sqrt(kappa*(weights@((g1['Pi']-g0['Pi'])**2)))/r).tolist(),
            W_RMS_gap=float(math.sqrt(weights@((g1['W']-g0['W'])**2))),
            F_RMS_gap=None if F0 is None or F1 is None else float(math.sqrt(weights@((F1-F0)**2)))))
    budget.check()
    return dict(common_points=len(fine.points),quadrature_order=list(order),
        metrics=[dict(metric_degree=L,classification=classification(e),**e.summary()) for L,e in zip(degrees,evaluations)],
        adjacent_expansion_and_field_gaps=gaps,physical_geometry_gaps=tensors,
        scope='The same finest-metric graph and points on matched observer-proper checkpoints and reference frames. Finite physical tensor/field samples, not a continuum or outermost-horizon bound.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records',type=Path,nargs='+',required=True)
    parser.add_argument('--cache-sha256',nargs='+')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--shape-degrees',type=int,nargs='+',default=[4,8])
    parser.add_argument('--snapshot',type=int,default=3)
    parser.add_argument('--budget-seconds',type=float,default=600.)
    parser.add_argument('--threads',type=int,choices=(1,2),default=2)
    parser.add_argument('--residual-tolerance-over-H',type=float,default=1e-6)
    parser.add_argument('--max-nfev',type=int,default=40)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    if not math.isfinite(args.budget_seconds) or args.budget_seconds<=0:raise ValueError('positive finite wall budget required')
    if not math.isfinite(args.residual_tolerance_over_H) or args.residual_tolerance_over_H<=0:raise ValueError('positive finite residual threshold required')
    if args.max_nfev<1 or args.snapshot<0:raise ValueError('nonnegative snapshot and positive optimizer budget required')
    shapes=args.shape_degrees
    if not shapes or any(d<0 for d in shapes) or any(b<=a for a,b in zip(shapes,shapes[1:])):
        raise ValueError('strictly increasing nonnegative shape degrees required')
    if args.cache_sha256 is not None and len(args.cache_sha256)!=len(args.records):
        raise ValueError('one caller cache pin is required per record, in matching order')
    pins=[None]*len(args.records) if args.cache_sha256 is None else args.cache_sha256
    sources=['scripts/check_nsc_atlas_valve_formation.py','src/recursive_horizons/nsc_atlas_trapped_parent.py',
        'src/recursive_horizons/nsc_atlas_trapped_surfaces.py']
    before={name:sha(ROOT/name) for name in sources}
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip();bound=True
    for name,expected in before.items():
        stored=subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True)
        bound=bound and stored.returncode==0 and hashlib.sha256(stored.stdout).hexdigest()==expected
    budget=Budget(args.budget_seconds)
    result=dict(schema='nsc-atlas-valve-formation-v1',status='in_progress',failure=None,
        budget_seconds=args.budget_seconds,threads=args.threads,selected_snapshot=args.snapshot,
        shape_degrees=shapes,initial_radii=list(START_RADII),residual_tolerance_over_H=args.residual_tolerance_over_H,
        producer_commit=commit,source_sha256=before,source_bound_to_commit=bound,
        evolution_runs=[],grid_diagnostics=[],surface_searches=[],common_surface_comparisons=[],future_trapping_candidates=[],
        physical_formation_admitted=False,universal_nonexistence_claim=False,
        admission_scope='Absent a producer cache pin, the explicitly supplied expected SHA is post-production caller admission, not preregistration or independent tamper-proof authentication.',
        focusing_scope='F=(Pi_S^2+Pi_chi^2+K_ijK^ij/kappa)/W at saved slices. It omits the normal acceleration-divergence term and is not itself a trapped-surface test.',
        search_scope='Finite star-shaped charts and starts on the finest supplied metric; both null signs and independent denser samples retained. Misses exclude no other surfaces, states or later times.')
    output.parent.mkdir(parents=True,exist_ok=True)
    def save():
        result['elapsed_seconds']=budget.elapsed
        result['producer_source_unchanged_during_run']=before=={name:sha(ROOT/name) for name in sources}
        output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    count=0
    def progress(message):
        nonlocal count
        count+=1;print(f'{count}. {message}',flush=True);save()
    def register_trapping(surface,evaluation,origin,center_name,metric_degree,initial_radius=None):
        if evaluation.future_trapped_sampled or evaluation.future_mots_sampled:
            row=candidate_summary((surface,evaluation,origin))
            row.update(center_name=center_name,center=surface.center.tolist(),metric_degree=metric_degree,
                snapshot=args.snapshot,initial_radius=initial_radius,physical_formation_admitted=False)
            result['future_trapping_candidates'].append(row)
    cohorts=[]
    try:
        set_num_threads(args.threads)
        with threadpool_limits(limits=args.threads),set_workers(args.threads):
            for path,pin in zip(args.records,pins):
                budget.check();path=(ROOT/path).resolve();cache=path.with_suffix('.npz')
                if not path.is_relative_to(ROOT) or not cache.is_relative_to(ROOT/'.build'):
                    raise ValueError('repository records with local regenerable state caches required')
                pin=None if pin is None else str(pin).lower()
                record=json.loads(path.read_bytes());producer_pin=record.get('state_cache_sha256')
                expected=producer_pin if producer_pin is not None else pin
                if expected is None or not re.fullmatch('[0-9a-f]{64}',str(expected)):
                    raise ValueError('an explicit caller SHA is required for each unpinned v2 cache')
                if pin is not None and pin!=expected:raise ValueError('caller and producer cache pins conflict')
                provider=GHSurfaceGeometry.from_state_cache(cache,record,expected_sha256=expected,snapshot=args.snapshot)
                with np.load(cache,allow_pickle=False) as data:times=data['times'].copy();states=data['states'].copy()
                selected={key:record[key] for key in ('degree','epsilon','status','stop_reason','prescribed_GH_gauge',
                    'gauge_constraint_damping','lapse_feedback_coefficient','temporal_tolerances','model_inputs','requested_observer_proper_endpoint') if key in record}
                selected.update(record_path=path.relative_to(ROOT).as_posix(),record_sha256=sha(path),
                    state_cache_path=cache.relative_to(ROOT).as_posix(),state_cache_sha256=expected,
                    cache_pin_origin='producer_record' if producer_pin is not None else 'post_production_caller_admission',
                    producer_commit=record['producer_commit'],source_sha256=record['source_sha256'],
                    selected_proper_duration=provider.proper_duration,selected_coordinate_t=provider.coordinate_t,
                    snapshot_constraints=[{key:point[key] for key in ('coordinate_t','observer_proper_duration',
                        'ADM_H_RMS_over_H_squared','ADM_M_RMS_over_H_squared','GH_C_RMS_over_H','metric_valid') if key in point}
                        for point in record['snapshots']])
                result['evolution_runs'].append(selected);cohorts.append((record,provider,times,states))
                progress(f'Authenticated metric degree {provider.degree}, {len(times)} saved events')
            cohorts.sort(key=lambda item:item[1].degree);degrees=[item[1].degree for item in cohorts]
            result['stabilization_threshold']=stabilization_threshold(cohorts[-1][1].system.parent)
            if len(set(degrees))!=len(degrees):raise ValueError('metric degrees must be distinct')
            common=cohorts[0][0]
            for record,provider,_,_ in cohorts:
                for key in ('epsilon','prescribed_GH_gauge','model_inputs'):
                    if record.get(key)!=common.get(key):raise ValueError('metric cohort has different physical preparation: '+key)
                if not math.isclose(provider.proper_duration,cohorts[0][1].proper_duration,rel_tol=0,abs_tol=2e-10):
                    raise ValueError('selected snapshots do not identify the same observer-proper event')
            for record,provider,times,states in cohorts:
                for index,(t,state,point) in enumerate(zip(times,states,record['snapshots'])):
                    row=grid_diagnostic(provider.system,state,t,index,point,budget);row['metric_degree']=provider.degree
                    result['grid_diagnostics'].append(row)
                    progress(f'Grid L={provider.degree}, event {index}, proper={row["observer_proper_duration"]:.6g}')
            providers=[LimitedGeometry(item[1],budget) for item in cohorts];fine=providers[-1];r=fine.r
            tolerance=args.residual_tolerance_over_H*r
            n=max(10,2*max(shapes)+4);validation=(2*n+1,4*n+3)
            result['common_validation_order']=list(validation)
            centers=(('scalar_anchor',anchor_point(fine.provider.system.parent)),('scalar_node',np.array([0.,0.,1.,0.])))
            for name,center in centers:
                candidates=[];search=dict(center_name=name,center=center.tolist(),metric_degree=degrees[-1],seeds=[],attempts=[])
                result['surface_searches'].append(search)
                for radius in START_RADII:
                    budget.check();seed=round_surface(center,radius,degree=shapes[0])
                    initial=evaluate_surface(seed,fine,frame_scale=r,quadrature_order=validation,expansion_tolerance=tolerance)
                    candidates.append((seed,initial,'round_seed'))
                    search['seeds'].append(dict(initial.summary(),initial_radius=radius))
                    register_trapping(seed,initial,'round_seed',name,degrees[-1],radius)
                    search['best_completed_candidate']=candidate_summary(min(candidates,key=lambda item:item[1].plus_max_abs))
                    progress(f'Start {name}, shape={shapes[0]}, radius={radius:g}')
                    found=find_mots(fine,center=center,initial_radius=radius,frame_scale=r,degree=shapes[0],
                        validation_order=validation,max_nfev=args.max_nfev,expansion_tolerance=tolerance)
                    summary=found.summary();summary.update(initial_radius=radius,classification=classification(found.evaluation))
                    search['attempts'].append(summary);candidates.append((found.surface,found.evaluation,'optimized'))
                    register_trapping(found.surface,found.evaluation,'optimized',name,degrees[-1],radius)
                    search['best_completed_candidate']=candidate_summary(min(candidates,key=lambda item:item[1].plus_max_abs))
                    progress(f'Root {name}, shape={shapes[0]}, max|theta+|/H={found.evaluation.plus_max_abs/r:.6g}, {summary["classification"]}')
                best=min(candidates,key=lambda item:item[1].plus_max_abs)
                for degree in shapes[1:]:
                    budget.check();warm=embed_degree(best[0],degree)
                    found=find_mots(fine,initial_surface=warm,frame_scale=r,validation_order=validation,
                        max_nfev=args.max_nfev,expansion_tolerance=tolerance)
                    summary=found.summary();summary.update(warm_start_shape_degree=best[0].degree,
                        classification=classification(found.evaluation))
                    search['attempts'].append(summary);candidates.append((found.surface,found.evaluation,'optimized'))
                    register_trapping(found.surface,found.evaluation,'optimized_warm_start',name,degrees[-1])
                    best=min(candidates,key=lambda item:item[1].plus_max_abs)
                    search['best_completed_candidate']=candidate_summary(best)
                    progress(f'Warm root {name}, shape={degree}, max|theta+|/H={found.evaluation.plus_max_abs/r:.6g}, {summary["classification"]}')
                surface,evaluation,origin=best
                search['selected_candidate']=candidate_summary(best)
                search['selected_candidate']['stress_diagnostic']=surface_diagnostic(evaluation,fine)
                comparison=compare_candidate(surface,providers,degrees,validation,tolerance,budget)
                comparison['center_name']=name;result['common_surface_comparisons'].append(comparison)
                for metric in comparison['metrics']:
                    if metric['future_trapped_sampled'] or metric['future_mots_sampled']:
                        row=dict(metric,center_name=name,center=surface.center.tolist(),snapshot=args.snapshot,
                            origin='common_finest_candidate',shape_degree=surface.degree,
                            shape_coefficients=surface.coefficients.tolist(),radius_bounds=list(surface.radius_bounds),
                            physical_formation_admitted=False)
                        result['future_trapping_candidates'].append(row)
                progress(f'Common {name} candidate evaluated on metric degrees {degrees}')
            result['status']='completed_finite_diagnostic'
    except TimeoutError as error:
        result['status']='budget_exhausted';result['failure']=str(error)
    except (ValueError,ArithmeticError,np.linalg.LinAlgError) as error:
        result['status']='failed_diagnostic';result['failure']=str(error)
    finally:
        if before!={name:sha(ROOT/name) for name in sources} and result['status']=='completed_finite_diagnostic':
            result['status']='source_changed_during_run';result['failure']='diagnostic source identity changed during execution'
        save()
    print(json.dumps(dict(status=result['status'],failure=result['failure'],elapsed_seconds=result['elapsed_seconds'],
        authenticated_runs=len(result['evolution_runs']),grid_events=len(result['grid_diagnostics']),
        searched_centers=len(result['surface_searches']),common_surface_comparisons=len(result['common_surface_comparisons']),
        universal_nonexistence_claim=False),sort_keys=True),flush=True)
    if result['status']!='completed_finite_diagnostic':
        raise SystemExit(2 if result['status']=='budget_exhausted' else 1)


if __name__=='__main__':main()
