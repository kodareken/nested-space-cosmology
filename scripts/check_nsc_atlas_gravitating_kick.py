"""Verify one finite kicked initial slice and its exact shear/coupling law."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_gravitating_kick import solve_kicked_initial_data,validate_kicked_data
from recursive_horizons.nsc_atlas_kick_relation import observer_relation

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree',type=int,default=16)
    parser.add_argument('--epsilon',type=float,default=.25)
    parser.add_argument('--kick-amplitude',type=float,default=1/32)
    parser.add_argument('--kick-source',type=int,choices=(1,2),default=1)
    parser.add_argument('--threads',type=int,default=4)
    parser.add_argument('--budget-seconds',type=float,default=180)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    names=['src/recursive_horizons/nsc_atlas_gravitating_kick.py','src/recursive_horizons/nsc_atlas_kick_york.py',
        'src/recursive_horizons/nsc_atlas_kick_relation.py','src/recursive_horizons/nsc_atlas_finite_parent_data.py',
        'src/recursive_horizons/nsc_atlas_s3_spectral.py','src/recursive_horizons/nsc_atlas_s3_frame.py',
        'src/recursive_horizons/nsc_atlas_parent_circulation.py','src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
        'scripts/check_nsc_atlas_gravitating_kick.py']
    def hashes():return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    before=hashes();commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bound=True
    for name in names:
        old=subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True)
        bound=bound and old.returncode==0 and hashlib.sha256(old.stdout).hexdigest()==before[name]
    started=time.perf_counter()
    with threadpool_limits(limits=args.threads):
        data=solve_kicked_initial_data(args.degree,epsilon=args.epsilon,kick_amplitude=args.kick_amplitude,
            source_harmonic=args.kick_source,budget_seconds=min(args.budget_seconds,120))
        validation=validate_kicked_data(data)
        relations=[observer_relation(epsilon=args.epsilon,source_harmonic=source,parent=data.parent) for source in (1,2)]
    elapsed=time.perf_counter()-started
    result=dict(schema='nsc-atlas-gravitating-kick-v1',
        status='finite_data_checked' if elapsed<=args.budget_seconds else 'budget_exhausted',
        epsilon=args.epsilon,kick_amplitude=args.kick_amplitude,kick_source=args.kick_source,
        initial_data=data.metadata,independent_constraint_validation=validation,observer_relations=relations,
        elapsed_seconds=elapsed,budget_seconds=args.budget_seconds,threads=args.threads,
        source_sha256=before,producer_commit=commit,source_bound_to_commit=bound,
        producer_source_unchanged_during_run=before==hashes(),
        exact_identity='Initially constant s, normal(s)=0 and trK=0 give proper d³ln(gKK)=-kappa mu² chi Pi_chi; the completed source has Pi=lambda f psi^-6 and Kij=lambda psi^-2 Abar.',
        scope='Exact momentum and conditional local relation, with a positive numerical Hamiltonian continuation. No unconditional nonlinear-branch theorem, physical TT-state selection, SI ruler or observational identification.')
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=result['status'],degree=args.degree,source=args.kick_source,
        validation=validation,ratios=[r['shear_to_logcoupling_jerk_ratio'] for r in relations],elapsed_seconds=elapsed),indent=2))
    if elapsed>args.budget_seconds:raise SystemExit(1)


if __name__=='__main__':main()
