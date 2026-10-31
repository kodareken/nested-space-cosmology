"""Actual-parent same-mode Dirac onset and bounded finite-slice controls."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.special import roots_jacobi
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_neck import solve_necks
from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve
from recursive_horizons.nsc_atlas_valve_dirac_modes import (
    initial_parent_mode_pump,finite_parent_mode_pump_check,radial_operator_blocks,full_mode_diagnostics,
)

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_valve_dirac_modes.py',
         'src/recursive_horizons/nsc_atlas_valve_spherical.py','src/recursive_horizons/nsc_atlas_valve_neck.py',
         'src/recursive_horizons/nsc_atlas_s3_spectral.py','src/recursive_horizons/nsc_atlas_parent_circulation.py',
         'src/recursive_horizons/nsc_atlas_parent_scale_selector.py','scripts/check_nsc_atlas_valve_dirac_onset.py')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def input_mode(initial):
    path=ROOT/'results/summaries/nsc-atlas-valve-dirac-modes-v1.json';record=json.loads(path.read_text())
    selected=next(v for v in record['levels'][-1]['lowest_by_sign'] if v['sign']>0)
    raw=ROOT/selected['raw_mode_cache']['path']
    if sha(raw)!=selected['raw_mode_cache']['sha256']:raise RuntimeError('input spatial-mode cache hash mismatch')
    with np.load(raw) as data:
        mode={key:data[key].copy() for key in ('degree','k','sign','eigenvalue','coefficients')}
        metadata=json.loads(str(data['metadata']))
    captured=record['provenance']['source_sha256_after']
    if metadata['source_sha256']!=captured:raise RuntimeError('input mode producer metadata mismatch')
    metric=ROOT/record['defined_geometry']['raw_cache']['path']
    if sha(metric)!=record['defined_geometry']['raw_cache']['sha256']:raise RuntimeError('input metric cache hash mismatch')
    with np.load(metric) as data:
        if not np.array_equal(initial.coefficients,data['coefficients']) or initial.parent.H!=float(data['H']):
            raise RuntimeError('initial geometry differs from the supplied spatial-mode geometry')
    archived=subprocess.check_output(['git','show','226d30e7:src/recursive_horizons/nsc_atlas_valve_neck.py'],cwd=ROOT)
    if hashlib.sha256(archived).hexdigest()!=captured['src/recursive_horizons/nsc_atlas_valve_neck.py']:
        raise RuntimeError('historical initial-neck source pin mismatch')
    # Re-evaluate the supplied mode, without selecting or fitting another mode.
    u,w=roots_jacobi(220,.5,.5);omega=initial.values(u)**2/initial.parent.H
    check=full_mode_diagnostics(*radial_operator_blocks(int(mode['degree']),int(mode['k']),u,w,omega),
                                u,w,omega,float(mode['eigenvalue']),mode['coefficients'])
    if check['full_PDE_residual_physical']>1e-6 or abs(check['physical_normalization']-1)>1e-12:
        raise RuntimeError('supplied spatial-mode operator/normalization check failed')
    return mode,float(selected['declared_mass']),dict(record=str(path.relative_to(ROOT)),record_sha256=sha(path),
        raw_mode=selected['raw_mode_cache'],defined_metric=record['defined_geometry']['raw_cache'],
        captured_input_source_pins=captured,historical_initial_neck_commit='226d30e7',
        current_operator_recheck=check,
        role='Captured input producer pins plus cache/geometry/operator checks; current appended onset analysis has separate source pins.')


def calculation(cache_directory):
    start=time.perf_counter();before={name:sha(ROOT/name) for name in SOURCES}
    with threadpool_limits(limits=1):
        initial=solve_necks(40,(1.,))[0];mode,mass,input_record=input_mode(initial)
        rows=[initial_parent_mode_pump(initial,mode,mass,orders) for orders in ((140,60),(220,100),(320,140))]
        system=SphericalValve(80);solution=system.evolve(.02)
        times=(.02,.01,.005,.0025);states=np.array([solution.sol(t).reshape(8,80) for t in times])
        finite=[]
        for t,state in zip(times,states):
            first=finite_parent_mode_pump_check(system,state,t,initial,mode,mass,isomap_order=128)
            second=finite_parent_mode_pump_check(system,state,t,initial,mode,mass,isomap_order=192)
            first['isomap_integral_order_refinement']=dict(orders=[128,192],
                lambda_value_spread=abs(first['lambda_value']-second['lambda_value']),
                lambda_ddot_spread=abs(first['lambda_ddot_finite']-second['lambda_ddot_finite']),
                pair_Hessian_spread=abs(first['pair_Hessian_finite']-second['pair_Hessian_finite']))
            first['pair_Hessian_difference_from_initial_limit']=first['pair_Hessian_finite']-rows[-1]['total_pair_Hessian']
            finite.append(first)
            print(t,first['pair_Hessian_finite'],first['pair_Hessian_difference_from_initial_limit'],flush=True)
        cache_directory=cache_directory.resolve();cache_directory.mkdir(parents=True,exist_ok=True)
        path=cache_directory/'finite-parent-slices.npz'
        np.savez_compressed(path,times=times,states=states,initial=system.initial,
                             metadata=json.dumps(dict(source_sha256=before,input_mode_record_sha256=input_record['record_sha256'],
                                                       preparation='Same80-mode one-mode parent; no source reset, through proper time0.02 only.')))
        massless=initial_parent_mode_pump(initial,mode,0.)
        round_initial=solve_necks(16,(0.,))[0];H=round_initial.parent.H
        round_mode=dict(degree=1,k=0,sign=1,eigenvalue=1.5*H,
                        coefficients=np.array([np.sqrt(H),0.,0.]))
        round_control=initial_parent_mode_pump(round_initial,round_mode,mass)
    after={name:sha(ROOT/name) for name in SOURCES}
    if before!=after:raise RuntimeError('source changed during bounded onset calculation')
    value=rows[-1]
    return dict(schema='nsc-atlas-valve-dirac-onset-v1',
                current_analysis=dict(source_sha256_before=before,source_sha256_after=after,source_unchanged=True,
                                      git_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),
                input_spatial_mode=input_record,initial_pump_quadrature=rows,
                quadrature_spreads=dict(total_pair_Hessian=abs(rows[-1]['total_pair_Hessian']-rows[-2]['total_pair_Hessian']),
                                       lambda_ddot=abs(rows[-1]['lambda_ddot']-rows[-2]['lambda_ddot']),
                                       projected_mass_ddot=abs(rows[-1]['projected_mass_ddot']-rows[-2]['projected_mass_ddot'])),
                finite_actual_parent_checks=finite,finite_slices_cache=dict(path=str(path.relative_to(ROOT)),sha256=sha(path)),
                independent_control=dict(massless_pair_Hessian=massless['total_pair_Hessian'],
                    round_de_sitter=dict(Omega_ddot_exact=H,lambda_ddot=round_control['lambda_ddot'],
                        lambda_ddot_expected=-H*H*round_mode['eigenvalue'],
                        projected_mass_ddot=round_control['projected_mass_ddot'],
                        pair_Hessian=round_control['total_pair_Hessian'],
                        pair_Hessian_expected=-mass*round_mode['eigenvalue']/round_control['energy']*H*H)),
                definitions=dict(fixed_initial_basis='beta_initial=-i*tau³*M_pn/6+O(tau4); coefficient is a propagator matrix element in the fixed initial canonical energy basis.',
                                 instantaneous_final_basis='beta_final=tau²*M_pn/(4E)+O(tau³) from final-basis rotation; not the cubic diagnostic.',
                                 singlet_normalization='For a normalized two-label singlet the fixed-basis probability coefficient isM_pn²/18, conditional on the low-mode initial vacuum.',
                                 coordinate_scope='Exact spatial slice isomap only. The time-dependent chart has a shift/unitary connection; only this same-lambda pair coefficient is shown unchanged.'),
                verdict=dict(nonzero_same_mode_initial_Hamiltonian_Hessian=value['total_pair_Hessian'],
                             fixed_initial_amplitude_cubic_coefficient=value['initial_reference_amplitude_cubic_coefficient'],
                             low_mode_initial_vacuum_is_an_added_postulate=True,
                             instantaneous_final_particle_number_computed=False,
                             UV_Hadamard_completion_selected=False,whole_field_occupation_selected=False,
                             renormalized_source_stress_computed=False,backreaction_computed=False,baryon_photon_eta=None,
                             scope='Connected intrinsic metric/radion/Dirac initial pump. Finite same-mode coefficient, not a full quantum state, detector observable or baryogenesis mechanism.'),
                numerical=dict(maximum_proper_time=.02,radial_order=220,isomap_orders=[128,192],spatial_modes=80,
                               finite_times=list(times),BLAS_threads=1,elapsed_seconds=time.perf_counter()-start,
                               caveat='Finite numerical and quadrature checks, not certified propagation of all initial metric/mode uncertainty; no higher-order or late-time extension.',
                               full_arrays_committed=False),
                reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_dirac_onset.py --output results/summaries/nsc-atlas-valve-dirac-onset-v1.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-dirac-onset-v1.json')
    parser.add_argument('--cache-directory',type=Path,default=ROOT/'.build/valve-dirac-onset-v1')
    args=parser.parse_args();result=calculation(args.cache_directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),initial=result['initial_pump_quadrature'][-1],
                          finite=result['finite_actual_parent_checks']),indent=2))


if __name__=='__main__':main()
