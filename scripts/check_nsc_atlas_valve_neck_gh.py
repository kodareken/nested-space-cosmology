#!/usr/bin/env python3
"""Independent full-GH control of the declared one-mode reflection neck.

The field preparation changes explicitly; the action and GH equations do
not. The initial state uses the package's production one-mode constraint
solver. Fresh run outputs default to .build; the previously measured
pre-promotion records are preserved separately and are not a rerun of this
byte-different producer.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from numba import set_num_threads
from scipy.fft import set_workers
from scipy.integrate import DOP853
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_gh import PAIRS,ParentGH
from recursive_horizons.nsc_atlas_s3_hopf import build_hopf_frame_basis
from recursive_horizons.nsc_atlas_trapped_parent import GHSurfaceGeometry
from recursive_horizons.nsc_atlas_trapped_surfaces import evaluate_surface,round_surface
from recursive_horizons.nsc_atlas_valve_neck import solve_necks,neck_diagnostics

ROOT=Path(__file__).resolve().parents[1]
VOLUME=2*math.pi**2
SOURCE_FILES=(
    'src/recursive_horizons/nsc_atlas_parent_gh.py',
    'src/recursive_horizons/nsc_atlas_s3_frame.py',
    'src/recursive_horizons/nsc_atlas_gh_pointwise.py',
    'src/recursive_horizons/nsc_atlas_gh_constraints.py',
    'src/recursive_horizons/nsc_atlas_finite_parent_data.py',
    'src/recursive_horizons/nsc_atlas_s3_hopf.py',
    'src/recursive_horizons/nsc_atlas_gh_observer.py',
    'src/recursive_horizons/nsc_atlas_parent_circulation.py',
    'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
    'src/recursive_horizons/nsc_atlas_parent_backreaction.py',
    'src/recursive_horizons/nsc_atlas_s3_spectral.py',
    'src/recursive_horizons/nsc_atlas_gh_matter_work.py',
    'src/recursive_horizons/nsc_atlas_passive_dop853.py',
    'src/recursive_horizons/nsc_atlas_gh_history.py',
    'src/recursive_horizons/nsc_atlas_kick_york.py',
    'src/recursive_horizons/nsc_atlas_gravitating_kick.py',
    'src/recursive_horizons/nsc_atlas_trapped_parent.py',
    'src/recursive_horizons/nsc_atlas_trapped_surfaces.py',
    'src/recursive_horizons/nsc_atlas_valve_neck.py',
)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def initial_state(system):
    """Production degree40 seed, projected into the unchanged full GH system."""
    z=system.parent.mass_squared/system.r**2
    normalization=24/(z+3)
    data=solve_necks(40,(1.,),parent=system.parent)[0]
    reference=neck_diagnostics(data)
    # Project psi^4, a polynomial of degree160, onto the full evolution basis.
    D=160+system.basis.degree
    startup=build_hopf_frame_basis(system.basis.degree,(math.ceil(D/4)+2,2*math.ceil((D+1)/2)))
    psi=data.values(startup.points[:,0])
    Fcoef=startup.project(psi**4)
    Fcoef[~np.all(system.basis.characters==[1,1],axis=1)]=0.
    F=startup.evaluate(Fcoef)
    if np.min(F)<=0:raise ValueError('projected one-mode initial metric lost positivity')
    C=np.zeros((24,len(system.basis.labels)))
    C[0,0]=-math.sqrt(VOLUME)
    for k,(a,b) in enumerate(PAIRS):
        if a==b and a:C[k]=Fcoef
    # Existing real l=1 basis has Y_(1,0,1,sin)=-sqrt(2)/pi X1.
    j=system.basis.labels.index((1,0,1,'sin'))
    amplitude=math.sqrt(normalization/system.parent.kappa)
    C[21,j]=-amplitude*math.pi/math.sqrt(2)
    for i,Di in enumerate(system.basis.generators):
        velocity=-.5*startup.evaluate(system.r*(Di@Fcoef))/F
        k=10+PAIRS.index((0,i+1));C[k]=startup.project(velocity)
        C[k,~np.all(system.basis.characters==system.chars[k],axis=1)]=0.
    # This observer sits on the exact X1=0 reflection neck. Its scalar
    # phase is identically undefined; only its geodesic proper age is used.
    system.X0=np.array([0.,1.,0.,0.]);system.update_point(system.X0)
    state=np.r_[system.pack(C),system.X0,[1.,0.,0.,0.],0.,0.]
    return state,dict(reference_initial_degree=40,reference_initial_row=reference,
        reference_initial_coefficients=data.coefficients.tolist(),
        projected_metric_degree=system.basis.degree,startup_quadrature_shape=list(startup.quadrature_shape),
        source='chi=A*X1; kappa*A^2=24/(z+3); all normal field momenta and K initially zero',
        A=amplitude,observer_initial_X=system.X0.tolist(),observer_clock='Initially normal geodesic on reflection neck; no scalar-phase clock.',
        metric_projection='Full psi40^4 projected, then exact finite-basis derivatives used for gauge velocity.',
        minimum_projected_F=float(np.min(F)))


def diagnostic(system,t,state):
    result=system.diagnostic(t,state)
    for name in ('observer_phase_wrapped','observer_sampled_unwrapped_phase'):
        result[name]=None
    provider=GHSurfaceGeometry(system,state,coordinate_t=t,
        provenance={'source_binding':'explicit_SHA256_bound_new_one_mode_preparation'})
    surface=round_surface([1.,0.,0.,0.],math.pi/2)
    cut=evaluate_surface(surface,provider,frame_scale=system.r,quadrature_order=(12,25))
    point=provider(state[system.dof:system.dof+4][None,:])
    result['neck_surface']=dict(cut.summary(),
        theta_plus_over_H_range=[float(np.min(cut.theta_plus)/system.r),float(np.max(cut.theta_plus)/system.r)],
        theta_minus_over_H_range=[float(np.min(cut.theta_minus)/system.r),float(np.max(cut.theta_minus)/system.r)],
        mean_curvature_max_abs_over_H=float(np.max(abs(cut.mean_curvature))/system.r),
        area_radius=float(math.sqrt(cut.area/(4*math.pi))),
        same_surface_expansion_difference_max_over_H=float(np.max(abs(cut.theta_plus-cut.theta_minus))/system.r),
        exact_reflection_surface='X1=0 fixed set, not a numerically searched radial zero')
    result['neck_normal_K_eigenvalues']=point['K_eigenvalues'][0].tolist()
    result['neck_lapse']=float(point['lapse'][0])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree',type=int,required=True)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--budget',type=float,default=180.)
    parser.add_argument('--proper',type=float,default=.02)
    parser.add_argument('--output',type=Path,
        help='Fresh path under .build; default uses the degree and a unique timestamp.')
    args=parser.parse_args()
    output=(ROOT/(args.output or Path(f'.build/valve-neck-gh-L{args.degree}-{time.time_ns()}.json'))).resolve()
    if output.exists() or not output.is_relative_to(ROOT/'.build'):raise ValueError('fresh output under .build required')
    output.parent.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    sources=SOURCE_FILES
    before={name:sha(ROOT/name) for name in sources}
    initializer={'script_sha256':sha(__file__),
        'one_mode_solver_sha256':sha(ROOT/'src/recursive_horizons/nsc_atlas_valve_neck.py')}
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    core_bound=True
    for name in sources:
        stored=subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True)
        core_bound=core_bound and stored.returncode==0 and hashlib.sha256(stored.stdout).hexdigest()==before[name]
    stored_script=subprocess.run(['git','show',commit+':scripts/check_nsc_atlas_valve_neck_gh.py'],cwd=ROOT,capture_output=True)
    initializer_bound=(stored_script.returncode==0 and hashlib.sha256(stored_script.stdout).hexdigest()==initializer['script_sha256']
        and core_bound)
    set_num_threads(args.threads)
    result=dict(schema='nsc-one-mode-neck-full-gh-control-v2',degree=args.degree,threads=args.threads,
        budget_seconds=args.budget,requested_neck_proper_duration=args.proper,producer_commit=commit,
        core_source_sha256=before,core_source_bound_to_commit=core_bound,
        initializer=initializer,initializer_bound_to_commit=initializer_bound,
        preparation_scope='New one-mode initial state, same existing action/normalization. Not the unchanged two-mode campaign.',
        physical_mode_suppression=False,SI_scale_selected=False,regular_child_proved=False)
    def guard():
        if time.perf_counter()-started>args.budget:raise TimeoutError('declared full GH control allowance exhausted')
    try:
        with threadpool_limits(limits=args.threads),set_workers(args.threads):
            system=ParentGH(args.degree,1.,gauge='lapse_controlled',backend='hopf')
            state,initial=initial_state(system);result['initial_preparation']=initial
            system.started=started;system.budget=args.budget
            guard();result['initial']=diagnostic(system,0.,state);guard()
            print('Initial L',args.degree,'Hres',result['initial']['ADM_H_RMS_over_H_squared'],
                  'Cres',result['initial']['GH_C_RMS_over_H'],'elapsed',time.perf_counter()-started,flush=True)
            def rhs(t,y):guard();return system.rhs(t,y)
            solver=DOP853(rhs,0.,state,t_bound=.1,rtol=3.e-9,atol=3.e-11,first_step=.002,max_step=.005)
            steps=0
            while solver.status=='running' and solver.y[system.dof+8]<args.proper:
                message=solver.step();steps+=1
                if solver.status=='failed':raise RuntimeError(message or 'DOP853 failure')
            if solver.y[system.dof+8]<args.proper:raise RuntimeError('proper endpoint not reached')
            dense=solver.dense_output()
            t=brentq(lambda value:dense(value)[system.dof+8]-args.proper,solver.t_old,solver.t,xtol=1.e-10)
            final=dense(t);guard();result['final']=diagnostic(system,t,final)
            result.update(status='endpoint_reached',steps=steps,rhs_calls=system.calls,
                proper_duration=float(final[system.dof+8]),coordinate_endpoint=t,
                neck_observer_X=final[system.dof:system.dof+4].tolist(),
                initial_to_final_observer_coordinate_displacement=float(np.linalg.norm(final[system.dof:system.dof+4]-system.X0)))
    except (TimeoutError,RuntimeError,ValueError,FloatingPointError) as exc:
        result.update(status='incomplete',stop_reason=str(exc))
    result['elapsed_seconds']=time.perf_counter()-started
    result['core_source_unchanged_during_run']=before=={name:sha(ROOT/name) for name in sources}
    result['initializer_unchanged_during_run']=initializer=={'script_sha256':sha(__file__),
        'one_mode_solver_sha256':sha(ROOT/'src/recursive_horizons/nsc_atlas_valve_neck.py')}
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(degree=args.degree,status=result['status'],elapsed=result['elapsed_seconds'],
        final=result.get('final'),output_sha256=sha(output)),sort_keys=True),flush=True)


if __name__=='__main__':main()
