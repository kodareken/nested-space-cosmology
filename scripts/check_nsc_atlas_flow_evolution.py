"""Bounded finite continuation and its leading fluid/gravity readouts."""
from concurrent.futures import ProcessPoolExecutor
import argparse,hashlib,json,time
from pathlib import Path
from recursive_horizons.nsc_atlas_flow_evolution import run_case
from recursive_horizons.nsc_atlas_flow_gravity import connections,boundary_coefficients
ROOT=Path(__file__).resolve().parents[1]


def worker(case):return run_case(*case,workers=4,cpu_limit=300.)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh repository output required')
    start=time.perf_counter();cases=[(64,.001,.05),(96,.001,.05),(128,.001,.05),(128,.0005,.05)]
    with ProcessPoolExecutor(max_workers=2) as pool:rows=list(pool.map(worker,cases))
    c,m,f,t=rows
    comparison={}
    for key in ('energy','dissipation'):
        vals=[row['final'][key] for row in rows]
        comparison[key]=dict(values=vals,coarse_change=abs(vals[1]-vals[0]),fine_change=abs(vals[2]-vals[1]),
            time_change=abs(vals[3]-vals[2]),empirical_error_indicator=2*max(abs(vals[2]-vals[1]),abs(vals[3]-vals[2])))
    comparison['generated_mean']=[]
    for i in range(2):
        vals=[row['generated_mean'][i] for row in rows]
        comparison['generated_mean'].append(dict(values=vals,coarse_change=abs(vals[1]-vals[0]),fine_change=abs(vals[2]-vals[1]),
            time_change=abs(vals[3]-vals[2]),empirical_error_indicator=2*max(abs(vals[2]-vals[1]),abs(vals[3]-vals[2]))))
    result=dict(schema='nsc-atlas-flow-evolution-v1',
        domain='Periodic cube [-2,2)^3, zero force, viscosity=.001, retained mode4 compact preparation, 0<=t<=.05.',
        algorithm='Strict 2/3-dealiased rotational Fourier Galerkin, orthogonal incompressible projection, SSPRK3; cubic periodic circle sampling.',
        cases=rows,refinement=comparison,
        gravity=connections(),
        final_boundary=boundary_coefficients(t['final']['energy'],t['final']['dissipation']),
        error_scope='Refinement/time-step differences are empirical numerical indicators, not rigorous continuum enclosures or a finite-amplitude bound on the gravity expansion.',
        energy_balance_max=max(abs(row['energy_balance_residual']) for row in rows),
        total_cpu_seconds=sum(row['cpu_seconds'] for row in rows),wall_seconds=time.perf_counter()-start,
        budget=dict(max_parallel_processes=2,FFT_workers_each=4,CPU_seconds_per_case=300,
                    cases_fixed_before_run=cases,no_grid_arrays_stored=True),
        physical_cosmology_or_baryon_fit=False)
    names=('src/recursive_horizons/nsc_atlas_flow_evolution.py','src/recursive_horizons/nsc_atlas_flow_gravity.py',
           'src/recursive_horizons/nsc_atlas_vortex_transfer.py','scripts/check_nsc_atlas_flow_evolution.py')
    result['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(refinement=comparison,energy_balance=result['energy_balance_max'],
                         CPU=result['total_cpu_seconds'],wall=result['wall_seconds']),indent=2))


if __name__=='__main__':main()
