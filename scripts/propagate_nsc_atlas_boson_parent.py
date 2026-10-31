"""Compute state-free coupled matter response on a resolved inherited parent."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from numba import set_num_threads
import scipy
from scipy.fft import set_workers
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_boson_hamiltonian import potential_hessian
from recursive_horizons.nsc_atlas_boson_parent import ParentBosonProbe
from recursive_horizons.nsc_atlas_boson_transport import propagate,symplectic_diagnostics
from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from record_nsc_atlas_gh_campaign import source_binding

ROOT=Path(__file__).resolve().parents[1]


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def response(probe,matrix,t):
    rows=probe.observer_rows(t);n=probe.n
    mass=np.diag(potential_hessian(np.zeros(2),sstar=probe.system.sstar))
    answers=[]
    for source in (1,2):
        z=matrix@probe.source_kick(source);q,p=z[:n],z[n:]
        norms=[]
        for field in (0,1):
            idx=probe.fields==field
            norms.append(float(math.sqrt(mass[field]*(q[idx]@q[idx])+p[idx]@p[idx])))
        answers.append(dict(source_harmonic=source,
            volume_reference_response_norm=norms[0],scalar_reference_response_norm=norms[1],
            volume_to_scalar_reference_response_ratio=norms[0]/norms[1] if norms[1]>0 else None,
            local_volume_field_per_unit_momentum_kick=float(rows['S']@q),
            local_scalar_field_per_unit_momentum_kick=float(rows['chi']@q),
            local_log_g_KK_per_unit_momentum_kick=float(rows['log_g_KK']@q)))
    return answers


def onset_response(probe):
    details=probe.hamiltonian_details(probe.initial_time)
    s=probe.s_indices;c=probe.chi_indices
    AS=details['A'][np.ix_(s,s)];AC=details['A'][np.ix_(c,c)]
    mixing=details['K'][np.ix_(s,c)]
    rows=probe.observer_rows(probe.initial_time)
    answer=[]
    for source in (1,2):
        f=probe.source_kick(source)[probe.n:][c]
        coefficient=-AS@mixing@AC@f/6
        answer.append(dict(source_harmonic=source,
            volume_q_cubic_coefficient_norm=float(np.linalg.norm(coefficient)),
            local_volume_cubic_coefficient=float(rows['S'][s]@coefficient),
            local_log_g_KK_cubic_coefficient=float(rows['log_g_KK'][s]@coefficient)))
    return answer


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history',type=Path,required=True)
    parser.add_argument('--background-record',type=Path,required=True)
    parser.add_argument('--probe-degree',type=int,required=True)
    parser.add_argument('--max-step',type=float,default=.02)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--budget-seconds',type=float,default=900)
    parser.add_argument('--decoupled-control',action='store_true')
    parser.add_argument('--matrix-cache',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    if not 1<=args.threads<=16 or not math.isfinite(args.budget_seconds) or args.budget_seconds<=0:raise ValueError('positive finite compute allocation required')
    history_path=(ROOT/args.history).resolve();background_path=(ROOT/args.background_record).resolve()
    if not history_path.is_relative_to(ROOT/'.build') or not background_path.is_relative_to(ROOT):
        raise ValueError('local regenerable history and repository background record required')
    cache=None
    if args.matrix_cache:
        cache=(ROOT/args.matrix_cache).resolve()
        if not cache.is_relative_to(ROOT/'.build') or cache.exists():raise ValueError('fresh matrix cache under .build required')
    background=json.loads(background_path.read_bytes())
    if background['status']!='endpoint_reached' or not source_binding(background):
        raise ValueError('verified source-bound background at its proper endpoint required')
    for name,expected in background['source_sha256'].items():
        if name.startswith('src/') and digest(ROOT/name)!=expected:
            raise ValueError('current model differs from the background producer: '+name)
    history_sha=digest(history_path)
    if history_sha!=background['dense_history']['sha256']:raise ValueError('dense history does not match its producer record')
    history=DenseGHHistory.load(history_path)
    if history.metadata['source_sha256']!=background['source_sha256'] or history.metadata['producer_commit']!=background['producer_commit']:
        raise ValueError('dense history source metadata mismatch')
    names=['src/recursive_horizons/nsc_atlas_boson_hamiltonian.py',
        'src/recursive_horizons/nsc_atlas_boson_parent.py','src/recursive_horizons/nsc_atlas_boson_transport.py',
        'src/recursive_horizons/nsc_atlas_gh_history.py','src/recursive_horizons/nsc_atlas_parent_gh.py',
        'src/recursive_horizons/nsc_atlas_s3_hopf.py','src/recursive_horizons/nsc_atlas_s3_frame.py',
        'scripts/propagate_nsc_atlas_boson_parent.py']
    names=list(dict.fromkeys(names+[name for name in background['source_sha256'] if name.startswith('src/')]))
    before={name:digest(ROOT/name) for name in names}
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bound=True
    for name in names:
        stored=subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True)
        bound=bound and stored.returncode==0 and hashlib.sha256(stored.stdout).hexdigest()==before[name]
    started=time.perf_counter();calls=0;result=None;failure=None;mass_min=math.inf
    set_num_threads(args.threads)
    with threadpool_limits(limits=args.threads),set_workers(args.threads):
        try:
            probe=ParentBosonProbe(history,args.probe_degree,decoupled_control=args.decoupled_control)
            def hamiltonian(t):
                nonlocal calls,mass_min
                if time.perf_counter()-started>args.budget_seconds:raise TimeoutError('bosonic transport compute budget exceeded')
                calls+=1
                if calls%100==0:print(f'Boson P={args.probe_degree} h={args.max_step:g}: t={t:.6g}/{probe.coordinate_endpoint:.6g}',flush=True)
                details=probe.hamiltonian_details(t)
                mass_min=min(mass_min,details['minimum_nodal_mass_eigenvalue'])
                return details['G']
            ends=[float(pair[1]) for pair in history.bounds if pair[1]<probe.coordinate_endpoint]
            edges=np.r_[probe.initial_time,ends,probe.coordinate_endpoint]
            result=propagate(hamiltonian,edges,max_step=args.max_step,checkpoints=probe.checkpoint_times)
            checkpoints=[]
            for t,matrix in zip(result.times,result.matrices):
                state=history.evaluate(float(t))
                checkpoints.append(dict(coordinate_t=float(t),
                    observer_proper_duration=float(state[probe.system.dof+8]),
                    symplectic=symplectic_diagnostics(matrix),response=response(probe,matrix,float(t))))
        except (TimeoutError,FloatingPointError,ValueError) as error:
            failure=str(error);result=None
    record=dict(schema='nsc-atlas-boson-response-v1',status='endpoint_reached' if result is not None else 'failed',
        failure=failure,probe_degree=args.probe_degree,max_step=args.max_step,threads=args.threads,
        budget_seconds=args.budget_seconds,elapsed_seconds=time.perf_counter()-started,
        hamiltonian_evaluations=calls,decoupled_control=args.decoupled_control,
        background_record_sha256=digest(background_path),background_history_sha256=history_sha,
        background_producer_commit=background['producer_commit'],background_degree=background['degree'],
        background_epsilon=background['epsilon'],background_temporal_tolerances=background['temporal_tolerances'],
        source_sha256=before,producer_commit=commit,source_bound_to_commit=bound,
        producer_source_unchanged_during_run=before=={name:digest(ROOT/name) for name in names},
        runtime=dict(python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__),
        inherited_covariance_selected=False,entropy_or_noise_selected=False,
        scope='Finite canonical Gaussian matter probe on the supplied nonlinear background; no metric fluctuations, quantum backreaction, continuum state or SI scale is selected.',
        response_normalization='Unit canonical momentum-density impulse in normalized source harmonic Y1 or Y2. Reference norms use the two vacuum masses; they are not physical energy-conversion fractions.',
        covariance_rule='Gamma(t)=M(t) Gamma_parent M(t)^T, retaining all cross correlations.',
        noise_rule='N(t,s)=P_q M(t) Gamma_parent M(s)^T P_q^T; incoming covariance remains an input.',
        retarded_rule='For Hamiltonian force -F q, G_R(t,s)=theta(t-s) P_q M(t) M(s)^(-1) P_p^T.')
    if result is not None:
        record.update(mode_counts=probe.counts,steps=result.steps,
            minimum_sampled_nodal_mass_eigenvalue=mass_min,
            onset_cubic_response=onset_response(probe),
            source_projection_residuals=probe.source_projection_residuals,checkpoints=checkpoints,
            max_checkpoint_symplectic_defect=result.max_checkpoint_symplectic_defect)
        if cache:
            np.savez_compressed(cache,times=result.times,matrices=result.matrices,fields=probe.fields,
                metadata=json.dumps(dict(probe_degree=args.probe_degree,history_sha256=history_sha,source_sha256=before)))
            record['matrix_cache_sha256']=digest(cache)
    output.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=record['status'],failure=failure,probe_degree=args.probe_degree,
        max_step=args.max_step,elapsed_seconds=record['elapsed_seconds'],
        final=record.get('checkpoints',[None])[-1]),indent=2))
    if result is None:raise SystemExit(1)


if __name__=='__main__':main()
