"""Reproduce actual initial-parent spatial Dirac modes and projector overlaps."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_neck import solve_necks,neck_diagnostics
from recursive_horizons.nsc_atlas_valve_source_gate import geometric_spinor_candidate
from recursive_horizons.nsc_atlas_valve_dirac_modes import calculate_dirac_level,omega_bounds

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_valve_dirac_modes.py',
         'src/recursive_horizons/nsc_atlas_valve_neck.py','src/recursive_horizons/nsc_atlas_s3_spectral.py',
         'src/recursive_horizons/nsc_atlas_parent_circulation.py','src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
         'src/recursive_horizons/nsc_atlas_inherited_charge.py','src/recursive_horizons/nsc_atlas_valve_source_gate.py',
         'scripts/check_nsc_atlas_valve_dirac_modes.py')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculate(cache_directory):
    start=time.perf_counter();before={name:sha(ROOT/name) for name in SOURCES}
    cache_directory=cache_directory.resolve();cache_directory.mkdir(parents=True,exist_ok=True)
    with threadpool_limits(limits=6):
        initial=solve_necks(40,(1.,))[0];candidate=geometric_spinor_candidate();mass=candidate['lowest_mass']
        metric=cache_directory/'defined-metric.npz'
        np.savez_compressed(metric,coefficients=initial.coefficients,degrees=np.arange(0,41,2),H=initial.parent.H,
                             epsilon=1.,degree=40,metadata=json.dumps(dict(source_sha256=before,
                                  definition='psi=sum(cl*C_l^1(u))/sqrt(2pi²); Omega=psi²/H')))
        metric_record=dict(path=str(metric.relative_to(ROOT)),sha256=sha(metric))
        levels=[]
        for degree in (8,16,24):
            level,modes=calculate_dirac_level(degree,initial=initial,mass=mass)
            for record,mode in zip(level['lowest_by_sign'],modes):
                path=cache_directory/f'degree-{degree}-sign-{mode["sign"]}.npz'
                np.savez_compressed(path,**mode,metadata=json.dumps(dict(source_sha256=before,
                                                                          defined_metric_sha256=metric_record['sha256'])))
                record['raw_mode_cache']=dict(path=str(path.relative_to(ROOT)),sha256=sha(path))
            levels.append(level)
            print(degree,[(r['sign'],r['eigenvalue'],r['full_PDE_residual_physical']) for r in level['lowest_by_sign']],flush=True)
        geometry_diagnostic=neck_diagnostics(initial)
    after={name:sha(ROOT/name) for name in SOURCES}
    if before!=after:raise RuntimeError('source changed during the bounded mode calculation')
    finest=levels[-1];positive=next(r for r in finest['lowest_by_sign'] if r['sign']>0)
    if max(v['full_PDE_residual_physical'] for v in finest['lowest_by_sign'])>1e-6:
        raise RuntimeError('fixed L24 full residual target was not reached')
    if not finest['compact_inverse_validation']['residual_interval_inside_isolated_global_first_cluster']:
        raise RuntimeError('the compact-inverse global-first isolation test failed')
    return dict(schema='nsc-atlas-valve-dirac-modes-v1',
                provenance=dict(source_sha256_before=before,source_sha256_after=after,source_unchanged=True,
                                git_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                                role='Current analysis source pins; raw coefficients are separately hashed, with no quantum state selected.'),
                defined_geometry=dict(rule='Same epsilon1 degree40 one-mode initial parent: h=Omega² gammaS3,Omega=psi40²/H.',
                                      isometry='Left quaternion multiplication Z=(-X4,X3,-X2,X1),Z4=oldX1; right frame preserved, not geometry averaging.',
                                      raw_cache=metric_record,Omega_bounds=omega_bounds(initial),
                                      Einstein_data_status='Defined finite psi polynomial; approximation to the Einstein constraint branch.',
                                      constraint_diagnostic=geometry_diagnostic),
                operator=dict(eigenproblem='D_gamma f=lambda Omega f; physical spinor=Omega^-1 f',
                              normalization='integral sqrt(1-u²)*Omega*(F²+G²)du=1; S2 spinor harmonics have unit angular norm.',
                              round_operator='D_gamma=i sigma_i D_i+3/2; existing right frame.',
                              sector='J=k+1/2 with multiplicity2k+2,orbital branchesj=k,k+1.',
                              exact_round_block='[[k+1.5,t],[t,-k-.5]],t=sqrt((l-k)(l+k+2)), after a constant i phase of branch2.',
                              residual='Direct positive integral of |D_gamma f-lambda Omega f|², also weighted byOmega^-1 for physical norm; includes output outside the trial space.',
                              residual_matrix_control='M2=P Omega² P is distinct from M@M. Subtraction can lose sign at roundoff; it is never used to fabricate a zero residual.'),
                declared_mass=dict(value=mass,source='Existing geometric_spinor_candidate lowest internal Dirac mass.',
                                   added_bulk_Dirac_probe=True,selected_by_spatial_geometry=False,
                                   energy_rule='E=sqrt(lambda²+mF0²); not an occupation or many-body source.'),
                levels=levels,
                projector=dict(definition='Normalized degeneracy-averaged first spatial projector; no occupied Fock state.',
                               density='rho_h=Omega^-2*(F²+G²)/(4pi)',
                               angular_control='For k0 and real F,G,f=(F I-iG sigma.n)xi/sqrt(4pi); f^dagger f=(F²+G²)I/(4pi).',
                               I4='integral |rho_h|²dV_h=(1/4pi)integral sqrt(1-u²)*Omega^-1*(F²+G²)²du',
                               actual_parent_volume=geometry_diagnostic['proper_total_volume'],
                               first_I4=positive['invariant_ground_projector_I4'],first_effective_volume=positive['effective_volume']),
                verdict=dict(first_cluster_identified_for_defined_finite_metric=True,
                             first_spatial_cluster_dimension=2,source_occupation_selected=False,
                             quantum_state_selected=False,contact_source_expectation_computed=False,
                             backreaction_computed=False,regular_child_proved=False,SI_scale_selected=False,
                             caveat='Analytic operator estimates with independently controlled numerical residuals; not interval-certified, nor a claim that the finite psi is the exact Einstein solution.'),
                numerical=dict(degrees=[8,16,24],BLAS_threads=6,elapsed_seconds=time.perf_counter()-start,
                               physical_mass_integrals='One-dimensional exact-polynomial Gauss-U rules; nonpolynomial residual/I4 independently checked at twice the radial order.',
                               full_arrays_committed=False),
                reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_dirac_modes.py --output results/summaries/nsc-atlas-valve-dirac-modes-v1.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-dirac-modes-v1.json')
    parser.add_argument('--cache-directory',type=Path,default=ROOT/'.build/valve-dirac-modes-v1')
    args=parser.parse_args();result=calculate(args.cache_directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),ground=result['levels'][-1]['lowest_by_sign'],
                          isolation=result['levels'][-1]['compact_inverse_validation']),indent=2))


if __name__=='__main__':main()
