"""Actual-parent initial multi-mode Dirac source operator; no state evolution."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import scipy
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_circulation import Parent
from recursive_horizons.nsc_atlas_valve_neck import NeckInitialData
import recursive_horizons.nsc_atlas_valve_dirac_modes  # Explicitly capture the operator's dynamic dependency.
from recursive_horizons.nsc_atlas_valve_source_gate import parent_dirac_source_operator,geometric_spinor_candidate

ROOT=Path(__file__).resolve().parents[1]
INPUTS=('results/summaries/nsc-atlas-valve-dirac-modes-v1.json',
        'results/summaries/nsc-atlas-valve-dirac-onset-v1.json')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def imported_project_sources():
    """Pin every loaded project Python source, including transitive owners."""
    package=(ROOT/'src/recursive_horizons').resolve();paths={Path(__file__).resolve()}
    for module in list(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if not name:continue
        path=Path(name).resolve()
        if path.suffix=='.py' and path.is_relative_to(package):paths.add(path)
    return {path.relative_to(ROOT).as_posix():sha(path) for path in sorted(paths)}


def load_inputs():
    pins={p:sha(ROOT/p) for p in INPUTS};spatial,onset=[json.loads((ROOT/p).read_text()) for p in INPUTS]
    metric=spatial['defined_geometry']['raw_cache'];path=ROOT/metric['path']
    if sha(path)!=metric['sha256']:raise RuntimeError('defined initial metric cache hash mismatch')
    with np.load(path) as raw:
        coeff=raw['coefficients'].copy();H=float(raw['H']);degree=int(raw['degree']);epsilon=float(raw['epsilon'])
    parent=Parent.retained()
    if parent.H!=H or degree!=40 or epsilon!=1.:raise RuntimeError('source operator requires the recorded one-mode psi40 epsilon1 parent')
    initial=NeckInitialData(parent,degree,epsilon,coeff,())
    for level in spatial['levels']:
        for mode in level['lowest_by_sign']:
            cache=mode['raw_mode_cache']
            if sha(ROOT/cache['path'])!=cache['sha256']:raise RuntimeError('preserved spatial mode cache hash mismatch')
    positive=next(v for v in spatial['levels'][-1]['lowest_by_sign'] if v['sign']>0)
    mass=float(positive['declared_mass'])
    if abs(mass-geometric_spinor_candidate()['lowest_mass'])>2e-14:raise RuntimeError('candidate mass differs from recorded initial mode')
    return initial,mass,pins,metric,onset['initial_pump_quadrature'][-1]


def calculation(cache_directory):
    start=time.perf_counter();before=imported_project_sources()
    initial,mass,input_pins,metric,onset=load_inputs()
    directory=cache_directory.resolve();directory.mkdir(parents=True,exist_ok=True)
    levels=[]
    with threadpool_limits(limits=6):
        for degree in (8,16,24):
            summary,matrices=parent_dirac_source_operator(initial,degree,mass=mass,k=0)
            path=directory/f'degree-{degree}-k-0.npz'
            np.savez_compressed(path,**matrices,metadata=json.dumps(dict(
                source_sha256=before,input_record_sha256=input_pins,defined_metric_sha256=metric['sha256'],
                degree=degree,k=0,mass=mass,preparation='Same defined initial psi40 epsilon1 parent.')))
            pair=matrices['positive_negative_pair_hessian'];index=int(matrices['ground_index'])
            summary['pair_matrix']=dict(shape=list(pair.shape),Frobenius_norm=float(np.linalg.norm(pair)),
                rows='Positive initial-energy labels with signed spatial eigenvalues.',
                columns='Negative initial-energy labels with signed spatial eigenvalues.',
                ground_column_index=index,
                ground_column_norm=math.sqrt(summary['ground_column_squared_coefficient_norm']),
                initial_reference_beta='-i*tau³*Hddot_pn/6+O(tau4), not an instantaneous-final particle spectrum.')
            summary['raw_matrices']=dict(path=str(path.relative_to(ROOT)),sha256=sha(path))
            levels.append(summary)
            print(degree,summary['same_mode_pair_Hessian'],summary['ground_column_squared_coefficient_norm'],
                  summary['other_mode_squared_coefficient_fraction'],flush=True)
        round_coeff=np.zeros_like(initial.coefficients);round_coeff[0]=math.sqrt(2*math.pi**2)
        round_initial=NeckInitialData(initial.parent,40,0.,round_coeff,())
        round_controls=[]
        for k in (0,1):
            control,matrices=parent_dirac_source_operator(round_initial,4,mass=mass,k=k)
            lam=control['ground_spatial_eigenvalue'];E=math.hypot(mass,lam)
            expected=-mass/E*initial.parent.H**2*lam
            pair=matrices['positive_negative_pair_hessian'];off=pair-np.diag(np.diag(pair))
            round_controls.append(dict(k=k,degree=4,same_mode_pair_Hessian=control['same_mode_pair_Hessian'],
                expected_same_mode_pair_Hessian=expected,absolute_difference=abs(control['same_mode_pair_Hessian']-expected),
                max_offdiagonal_pair_Hessian=float(np.max(abs(off))),
                Gram_identity_error=control['Gram_identity_error'],
                original_Dirac_projection_error=control['original_Dirac_projection_error'],
                other_mode_squared_coefficient_fraction=control['other_mode_squared_coefficient_fraction'],
                role='Round deSitter fixed-operator control only.'))
        massless,matrices=parent_dirac_source_operator(initial,4,mass=0.,k=0)
        massless_control=dict(degree=4,k=0,same_mode_pair_Hessian=massless['same_mode_pair_Hessian'],
            ground_column_squared_coefficient_norm=massless['ground_column_squared_coefficient_norm'],
            max_absolute_pair_Hessian=float(np.max(abs(matrices['positive_negative_pair_hessian']))),
            physical_massless_candidate_selected=False,
            role='Hypothetical external mass0 operator control: same-mode zero need not imply intermode matrix zero. The untwisted internal candidate remains massive.')
    after=imported_project_sources()
    if before!=after or input_pins!={p:sha(ROOT/p) for p in INPUTS}:
        raise RuntimeError('actual imported project sources or preserved input records changed during calculation')
    fine=levels[-1]
    comparisons={key:dict(values=[v[key] for v in levels],adjacent_spreads=[abs(levels[i+1][key]-levels[i][key]) for i in (0,1)])
                 for key in ('same_mode_pair_Hessian','ground_column_squared_coefficient_norm','other_mode_squared_coefficient_fraction')}
    return dict(schema='nsc-atlas-valve-dirac-source-v1',
                provenance=dict(source_sha256_before=before,source_sha256_after=after,source_unchanged=True,
                    input_record_sha256=input_pins,defined_metric_cache=metric,
                    git_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                    role='Current initial source-operator analysis; older input producer pins and data are preserved separately.'),
                preparation=dict(initial_degree=40,epsilon=1.,mass=mass,k=0,
                    geometry='The same defined finite one-mode initial parent; approximation to the Einstein constraint branch.',
                    domain='Initial proper chart only. No R grid or time evolution.'),
                operator=dict(meaning='Full finite-k initial radial Dirac Hamiltonian Hessian with all spatial modes in its declared cutoff.',
                    units='Derivative of Hhat=Lunit*Hphysical with respect to tauhat=tauphysical/Lunit, so Mhat=Lunit^3*Mphysical; squared coefficient norms use these retained dimensionless units. No SI ruler is selected.',
                    chart='Fixed original proper chart and canonical half density; no time-dependent isotropic-chart shift omitted.',
                    eigenbasis='Both signed spatial Dirac eigenvalues in the k0 block; positive/negative INITIAL energy spinors at each label.',
                    checks='Polar Gram versus Omega-weighted norm, reconstructed original Dirac matrix, Hessian hermiticity and independent same-mode isotropic Hellmann-Feynman value.',
                    fractions='Squared matrix-amplitude-coefficient column norms only, not energy fractions, probabilities, occupation data or a complete UV source.'),
                levels=levels,adjacent_spreads=comparisons,
                independent_prior_onset=dict(record_value=onset['total_pair_Hessian'],
                    finest_operator_value=fine['same_mode_pair_Hessian'],
                    absolute_difference=abs(fine['same_mode_pair_Hessian']-onset['total_pair_Hessian'])),
                round_de_sitter_controls=round_controls,massless_operator_control=massless_control,
                verdict=dict(same_mode_onset_recovered=abs(fine['same_mode_pair_Hessian']-onset['total_pair_Hessian'])<1e-10,
                    nonzero_other_mode_operator_entries=fine['other_mode_squared_coefficient_fraction']>1e-10,
                    same_mode_closure_is_not_the_full_quadratic_operator=fine['other_mode_squared_coefficient_fraction']>1e-10,
                    whole_field_initial_state_selected=False,UV_Hadamard_completion_selected=False,
                    renormalized_stress_computed=False,interaction_connected_four_point_computed=False,
                    g_selected=False,backreaction_computed=False,baryon_photon_eta=None,
                    scope='New operator/preparation connection, not a physical particle spectrum, full kinetic stress or turnaround prediction.'),
                numerical=dict(degrees=[8,16,24],k=0,BLAS_threads=6,
                    numpy_version=np.__version__,scipy_version=scipy.__version__,
                    elapsed_seconds=time.perf_counter()-start,full_arrays_committed=False),
                reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_dirac_source.py --output results/summaries/nsc-atlas-valve-dirac-source-v1.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-dirac-source-v1.json')
    parser.add_argument('--cache-directory',type=Path,default=ROOT/'.build/valve-dirac-source-v1')
    args=parser.parse_args();result=calculation(args.cache_directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),finest=result['levels'][-1],
                          round=result['round_de_sitter_controls'],massless=result['massless_operator_control']),indent=2))


if __name__=='__main__':main()
