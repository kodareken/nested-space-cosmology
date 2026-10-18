"""Same-field wall/ray precision diagnosis after the failed mass derivative."""
from __future__ import annotations

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['NUMBA_NUM_THREADS']='1'

import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_compact_source import CompactMassChart
from recursive_horizons.nsc_atlas_dirac_mean import PreparedDiracSource
from recursive_horizons.nsc_atlas_dirac_backreaction import MeanBackreaction
from measure_nsc_atlas_dirac_mass_response import SOURCES

ROOT=Path(__file__).resolve().parents[1]


def diagnose(mass):
    with threadpool_limits(limits=1):
        started=time.process_time()
        chart=CompactMassChart(WallRadarChart(),mass)
        source=PreparedDiracSource(chart,128)
        U,V=.8*chart.duration,.7*chart.duration
        def remaining():
            value=50-(time.process_time()-started)
            if value<=0:
                raise RuntimeError('Same-field diagnosis exhausted 50 CPU seconds')
            return value
        coarse=MeanBackreaction(source,U,wall_cpu_limit=remaining())
        tight=MeanBackreaction(source,U,wall_rtol=2e-11,wall_atol=2e-13,
            wall_cpu_limit=remaining())
        variants=[]
        for name,model,rtol,atol in (
                ('original',coarse,2e-8,2e-10),
                ('tight_wall_only',tight,2e-8,2e-10),
                ('tight_wall_and_ray',tight,2e-11,2e-13)):
            row=model.clock_response(U,V,order=12,ray_rtol=rtol,ray_atol=atol,
                cpu_limit=remaining())
            wall=model.wall_state(U)
            # Independent direct integral of the two wall-mass source laws;
            # the ODE interpolants are not used as the mass integrand.
            def derivative(t):
                data=source.evaluate(t,t)['orthonormal']
                g=chart.wall_geometry(t)
                j,pr=data['outward_energy_flux'],data['radial_pressure']
                R=g['x']
                return -.5*R*R*np.array([g['beta_plus']*j+g['v']*pr,
                                        g['beta_minus']*j+g['v']*pr])
            direct=model._integrate(source.start,U,derivative,12)
            variants.append(dict(variant=name,clock=row['frequency_fractional_change'],
                wall_exterior_mass=wall['exterior_mass'],wall_cap_mass=wall['cap_mass'],
                wall_displacement=wall['displacement'],
                direct_wall_mass_gaps=[wall['exterior_mass']-float(direct[0]),
                                       wall['cap_mass']-float(direct[1])],
                normalized_bulk_mass_path_gap=abs(row['mass_path_gap'])/max(1.,abs(row['mass'])),
                wall_rhs=model.wall_nfev,ray_rhs=row['ray_rhs_calls']))
        return dict(central_mass=mass,cells=128,variants=variants,pid=os.getpid(),
            process_CPU_total=time.process_time())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('use an unused compact diagnosis path')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    hashes={}
    for name in (*SOURCES,'scripts/diagnose_nsc_atlas_dirac_source_precision.py'):
        raw=(ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT)!=raw:
            raise ValueError('freeze precision diagnosis before execution: '+name)
        hashes[name]=hashlib.sha256(raw).hexdigest()
    started=time.perf_counter()
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        jobs=[pool.submit(diagnose,b) for b in (0.,5e-5)]
        rows=[future.result() for future in as_completed(jobs)]
    rows.sort(key=lambda r:r['central_mass'])
    by_name={}
    for i in range(3):
        name=rows[0]['variants'][i]['variant']
        by_name[name]=(rows[1]['variants'][i]['clock']-rows[0]['variants'][i]['clock'])/5e-5
    charged=time.process_time()+sum(r['process_CPU_total'] for r in rows)+4
    record=dict(schema='nsc-atlas-dirac-source-precision-v1',source_commit=commit,
        source_sha256=hashes,rows=rows,one_sided_mass_derivatives=by_name,
        runtime=dict(wall_seconds=time.perf_counter()-started,charged_CPU_seconds=charged,
            process_workers=2,per_worker_threads=1,new_field_resolutions=0,new_P001_grids=0),
        limits=dict(charged_CPU_seconds=120,summary_bytes=12288),
        scope='Two existing N=128 fields reused across original/tight wall and ray tolerances; '
              'direct wall-mass quadrature diagnoses the failed source derivative. No new source or grids.')
    if charged>120:
        record['status']='precision_diagnosis_budget_exhausted'
    else:
        record['status']='same_field_precision_diagnosis_evaluated'
    raw=json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>12288:
        raise RuntimeError('Precision diagnosis exceeds compact summary allowance')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps(record),flush=True)
    if charged>120:
        raise SystemExit(1)


if __name__=='__main__':
    main()
