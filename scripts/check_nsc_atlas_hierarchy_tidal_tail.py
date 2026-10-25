"""Diagnose the exterior quadrupole and test its harmonic phase boundary."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from recursive_horizons.nsc_atlas_hierarchy_perturbations import PerturbationRun
from recursive_horizons.nsc_atlas_hierarchy_tidal_tail import run_tidal_tail
from check_nsc_atlas_hierarchy_perturbations import INPUTS as BULK_INPUTS

ROOT=Path(__file__).resolve().parents[1]
INPUTS=BULK_INPUTS+('src/recursive_horizons/nsc_atlas_hierarchy_tidal_tail.py',)


def worker(job):
    name,config=job
    key=hashlib.sha256(json.dumps(dict(config=asdict(config),sources={
        f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}),sort_keys=True).encode()).hexdigest()
    folder=ROOT/'.build/hierarchy-tidal-tail'/key;path=folder/'completed.json'
    if path.is_file():r=json.loads(path.read_text(encoding='utf-8'))
    else:
        folder.mkdir(parents=True,exist_ok=True);r,raw=run_tidal_tail(config,retain=True)
        np.savez_compressed(folder/'fields.npz',radius=raw['radius'],final_u=raw['final_u'],final_q=raw['final_q'])
        (folder/'history.json').write_text(json.dumps(raw['history']),encoding='utf-8')
        path.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    return name,r


def calculation():
    wall=time.perf_counter();cpu=time.process_time();base=PerturbationRun()
    configs=dict(coarse=replace(base,points=4095,step=.04),medium=base,
        outer_control=replace(base,outer=1536.,points=12287))
    with ProcessPoolExecutor(max_workers=3) as pool:runs=dict(pool.map(worker,configs.items()))
    old=json.loads((ROOT/'results/summaries/nsc-atlas-hierarchy-perturbations-v1.json').read_text(encoding='utf-8'))
    source={name:r['source_final_field_sha256']==old['runs'][name]['source_final_field_sha256'] for name,r in runs.items()}
    key='potential_RMS_coefficient_per_initial_fractional_field_norm'
    value=runs['medium'][key]
    grid=abs(value-runs['coarse'][key]);domain=abs(value-runs['outer_control'][key])
    error=2*max(grid,domain)
    errors=[]
    if not all(source.values()):errors.append('spherical source changed')
    if max(r['maximum_symplectic_gap'] for r in runs.values())>1e-8:errors.append('Hamiltonian boundary transport failed')
    total=time.process_time()-cpu+sum(r['cpu_seconds'] for r in runs.values())
    result=dict(schema='nsc-atlas-hierarchy-tidal-tail-v1',
        result='Leading homogeneous-ambient tidal-phase matching for the external quadrupole, following an unresolved Dirichlet probe.',
        runs=runs,source_replay_exact_matches=source,
        final_potential_RMS_coefficient=value,numerical_comparison_envelope=error,
        grid_change=grid,outer_domain_change=domain,
        outward_probe_resolved=bool(error<.1*abs(value)),
        old_Dirichlet_probe=old['external_quadrupole_potential'],
        limits='Separate leading harmonic-tail boundary model. The no-incoming-density harmonic exterior has eta proportional to i r^-3 and q proportional to i r^-2. General outgoing radiative and Jeans responses are not solved by this local condition. A finite-domain comparison tests its usefulness here.',
        numerical_validation_errors=errors,process_cpu_seconds_including_workers=total,
        compute_budget_seconds=800,budget_exhausted=total>800,wall_seconds=time.perf_counter()-wall,
        timing_scope='CPU figures include the underlying computations reused from authenticated checkpoints; wall_seconds is this invocation.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS+('scripts/check_nsc_atlas_hierarchy_tidal_tail.py',
            'tests/test_nsc_atlas_hierarchy_tidal_tail.py','results/summaries/nsc-atlas-hierarchy-perturbations-v1.json')},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_hierarchy_tidal_tail.py --output results/summaries/nsc-atlas-hierarchy-tidal-tail-v1.json')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);args=p.parse_args();r=calculation()
    encoded=json.dumps(r,indent=2,allow_nan=False)+'\n'
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded)
    print(json.dumps(dict(value=r['final_potential_RMS_coefficient'],error=r['numerical_comparison_envelope'],
        domain_change=r['outer_domain_change'],resolved=r['outward_probe_resolved'],cpu=r['process_cpu_seconds_including_workers'])))
    if r['numerical_validation_errors']:raise RuntimeError(r['numerical_validation_errors'])
    if r['budget_exhausted']:raise RuntimeError('800 CPU-second tail-matching budget exhausted; diagnose before further computation')
