"""Reproduce the same-state local conformal continuation, overlap and controls."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_spherical import (
    SphericalValve,ConformalValve,align_conformal_center,
    conformal_de_sitter_control,conformal_forcing,
)


ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_valve_spherical.py',
         'src/recursive_horizons/nsc_atlas_valve_neck.py',
         'src/recursive_horizons/nsc_atlas_parent_circulation.py',
         'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
         'scripts/check_nsc_atlas_valve_conformal.py')
ENDPOINT=1.1481865519931276


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_cache(path,**arrays):
    path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,**arrays)
    return dict(path=str(path.relative_to(ROOT)),sha256=sha(path))


def comparison(values):
    differences=[abs(values[i+1]-values[i]) for i in (0,1)]
    return dict(values=values,adjacent_differences=differences,
                difference_ratio=differences[0]/differences[1] if differences[1] else None,
                empirical_fine_grid_error_indicator=differences[1],
                meaning='Full successive-grid difference, not a certified continuum bound.')


def calculation(cache_directory):
    start=time.perf_counter();cache_directory=cache_directory.resolve()
    before={name:sha(ROOT/name) for name in SOURCES}
    git_parent=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    # The archived original module belongs to that prior scientific version.
    # The implementation being evaluated below belongs to current working-tree pins.
    archived=subprocess.check_output(['git','show','8dbe971d:src/recursive_horizons/nsc_atlas_valve_spherical.py'],cwd=ROOT)
    with threadpool_limits(limits=1):
        system=SphericalValve(80,initial_degree=40,budget_seconds=60.)
        solution=system.evolve(2.);state=solution.y[:,-1].reshape(8,80)
        handoff_diagnostic=system.diagnostic(2.,state)
        handoff_cache=save_cache(cache_directory/'handoff.npz',proper_time=2.,state=state,u=system.u,
                                 metadata=json.dumps(dict(source_sha256=before,proper_time=2.,modes=80,initial_degree=40)))
        frame=ConformalValve(system,state,initial_proper_time=2.,u_max=.9)
        # This is an overlap-only continuation of the same reference state.
        # It neither replaces nor extends the later conformal calculation.
        reference=solve_ivp(system.rhs,(2.,2.2),state.ravel(),method='DOP853',
                            rtol=3e-11,atol=1e-13,max_step=.002)
        if not reference.success:raise RuntimeError(reference.message)
        reference_state=reference.y[:,-1].reshape(8,80)
        reference_diagnostic=system.diagnostic(2.2,reference_state)
        reference_fields=system.pieces(reference_state,[0.])
        reference_values=dict(center_radius=float(reference_fields['radius'][0]),
                              center_S=float(reference_fields['y'][4,0]),
                              center_Pi_S=float(reference_fields['y'][6,0]),
                              center_theta_plus_over_H=-2*float(reference_fields['kb'][0])/system.H)
        levels=[];results=[];alignments=[]
        for points in (2001,4001,8001):
            result=frame.evolve(points,ENDPOINT,budget_seconds=60.)
            if result['status']!='requested_endpoint':raise RuntimeError((points,result['status'],result['message']))
            results.append(result);profile=result['end_profile']
            raw=save_cache(cache_directory/f'points-{points}.npz',
                           center_history=result['center_history'],center_columns=np.asarray(result['center_columns']),
                           **profile,metadata=json.dumps(dict(source_sha256=before,points=points,
                                                               handoff_cache_sha256=handoff_cache['sha256'])))
            levels.append(dict(points=points,step=result['step'],status=result['status'],
                               initial=result['rows'][0],endpoint=result['rows'][-1],
                               events=result['events'],raw_cache=raw,
                               maximum_implicit_iterations=result['maximum_implicit_iterations'],
                               elapsed_seconds=result['elapsed_seconds']))
            alignments.append(align_conformal_center(result,2.2))
            print(points,'tau',profile['center_proper_time'],'R',result['rows'][-1]['center_radius'],
                  'Ham/H2',result['rows'][-1]['Hamiltonian_RMS_over_H2'],flush=True)
        endpoint_rows=[v['endpoint'] for v in levels]
        if max(r['conformal_time'] for r in endpoint_rows)-min(r['conformal_time'] for r in endpoint_rows)>1e-10:
            raise RuntimeError('the spatial refinements did not reach one common conformal event')
        fields=('center_proper_time','center_radius','center_Rproper','center_Rproper_acc',
                'center_theta_plus_over_H','center_theta_product_over_H2','center_R_gradient_norm_squared',
                'center_spacetime_Hess_ss_over_R_H2','center_Rll_over_R_H2','Hamiltonian_RMS_over_H2')
        refinements={field:comparison([r[field] for r in endpoint_rows]) for field in fields}
        overlap={field:dict(synchronous_value=value,
                            conformal_values=[v['values'][field] for v in alignments],
                            interpolation_indicators=[v['interpolation_indicator'][field] for v in alignments],
                            finest_difference=alignments[-1]['values'][field]-value)
                 for field,value in reference_values.items()}
        controls=[conformal_de_sitter_control(n,system=system) for n in (501,1001,2001)]
        H=system.H;T=.3;x=np.linspace(-.7,.7,101);scale=1/math.cos(H*T)
        R=scale*np.cos(H*x)/H;a=-math.log(math.cos(H*T));aT=H*math.tan(H*T);zero=np.zeros_like(x)
        y=np.array([R,np.full_like(x,a),zero,zero]);yt=np.array([aT*R,np.full_like(x,aT),zero,zero])
        yx=np.array([-scale*np.sin(H*x),zero,zero,zero])
        exact=np.array([(H*H*scale*scale+aT*aT+H*H)*R,np.full_like(x,H*H*scale*scale),zero,zero])
        identity=float(np.max(abs(conformal_forcing(y,yt,yx,system.sstar,system.kappa)-exact)))
    after={name:sha(ROOT/name) for name in SOURCES}
    if before!=after:raise RuntimeError('source changed during the calculation; do not publish mixed-source results')
    return dict(schema='nsc-atlas-valve-conformal-v1',
                provenance=dict(role='Current working-tree source pins; not the archived prior Git version.',
                                source_sha256_before=before,source_sha256_after=after,source_unchanged=True,
                                git_parent_commit=git_parent,archived_prior_parent_commit='8dbe971d',
                                archived_prior_module_sha256=hashlib.sha256(archived).hexdigest()),
                handoff=dict(proper_time=2.,modes=80,initial_constraint_degree=40,epsilon=1.,
                             diagnostic=handoff_diagnostic,raw_cache=handoff_cache,
                             initial_radial_half_width=frame.half_width,initial_u_max=.9,
                             rule='a=0,a_T=-K_rad,R_T=-K_tan R,phi_T=Pi; X is proper radial distance. Fields and momenta are retained.'),
                equations=dict(metric='e^(2a)(-dT²+dX²)+R²dOmega²',
                               R_wave='R_TT-R_XX=-(R_T²-R_X²)/R-e^(2a)/R+kappa*W*R*e^(2a)',
                               a_wave='a_TT-a_XX=(R_T²-R_X²)/R²+e^(2a)/R²-kappa*sum(phi_T²-phi_X²)/2',
                               field_wave='phi_TT-phi_XX=-2*(R_T phi_T-R_X phi_X)/R-e^(2a)*W_phi',
                               boundary='Shrinking characteristic Cauchy diamond; no new boundary state.',
                               Hessian='Hess_ss R=R_ll+K_rad R_proper; R_ll alone depends on slicing.'),
                common_conformal_endpoint=endpoint_rows[-1]['conformal_time'],levels=levels,
                refinement=refinements,
                overlap=dict(proper_time=2.2,synchronous_diagnostic=reference_diagnostic,observables=overlap,
                             clock='Center proper time integral; same reflection-center geodesic.',
                             interpolation='Actual-step PCHIP; full-vs-decimated history change reported as an empirical indicator.'),
                exact_de_sitter=dict(levels=controls,analytic_PDE_identity_max_abs_residual=identity,
                                     radius_difference_ratios=[controls[i]['max_radius_error']/controls[i+1]['max_radius_error'] for i in (0,1)],
                                     metric_difference_ratios=[controls[i]['max_metric_log_error']/controls[i+1]['max_metric_log_error'] for i in (0,1)],
                                     constraint_noise='Small finite-difference constraint controls are reported honestly, including nonmonotone noise; field convergence is the analytic-control comparison.'),
                numerical=dict(method='Second-order implicit-centered characteristic diamond; CFL1; initial Cauchy Taylor order2',
                               proper_clock='Trapezoidal lapse integral',BLAS_threads=1,
                               requested_conformal_endpoint=ENDPOINT,elapsed_seconds=time.perf_counter()-start,
                               full_arrays_committed=False),
                verdict=dict(same_state_continues_beyond_old_global_numerical_break=True,
                             angular_contraction_persists=True,temporal_bounce_observed=False,
                             all_retained_patch_spheres_future_trapped=endpoint_rows[-1]['full_current_patch_future_trapped'],
                             regular_child_proved=False,physical_curvature_singularity_proved=False,
                             source_state_selected_by_inheritance=False,
                             scope='Finite empirical numerical errors, not a continuum enclosure. The local causal edge is not the full receiving-band marginal boundary; late whole-band influx is not inferred. Zero center flux permits symmetric feeding.',
                             next='Keep the finite continuation fixed; derive the next interior connection in the same action.'),
                reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_conformal.py --output results/summaries/nsc-atlas-valve-conformal-v1.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-conformal-v1.json')
    parser.add_argument('--cache-directory',type=Path,default=ROOT/'.build/valve-conformal-v1')
    args=parser.parse_args();result=calculation(args.cache_directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),endpoint=result['levels'][-1]['endpoint'],
                          radius_refinement=result['refinement']['center_radius'],overlap=result['overlap']['observables']),indent=2))


if __name__=='__main__':main()
