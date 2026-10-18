"""Bounded clock-energy state candidate on the actual nested Dirac geometry.

Finite radial/angular cutoffs are explicit. This is neither a renormalized
vacuum stress nor a full scalar-potential slip calculation.
"""
from __future__ import annotations

import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
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
from scipy import constants
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_compact_source import CompactMassChart
from recursive_horizons.nsc_atlas_dirac_preparation import ClockEnergyBasis

ROOT=Path(__file__).resolve().parents[1]
CPU_LIMIT=240.
SOURCES=(
    'src/recursive_horizons/nsc_atlas_dirac_preparation.py',
    'src/recursive_horizons/nsc_atlas_dirac_mean.py',
    'src/recursive_horizons/nsc_atlas_dirac_stress.py',
    'src/recursive_horizons/nsc_atlas_dirac_radar.py',
    'src/recursive_horizons/nsc_atlas_dirac_compact_source.py',
    'src/recursive_horizons/nsc_atlas_dirac_scattering.py',
    'src/recursive_horizons/nsc_atlas_cap_wall_response.py',
    'src/recursive_horizons/nsc_atlas_causal_reference.py',
    'src/recursive_horizons/nsc_atlas_dirac_evolution.py',
    'src/recursive_horizons/nsc_influence.py',
    'tests/test_nsc_atlas_dirac_preparation.py',
    'tests/test_nsc_atlas_dirac_mean.py',
    'scripts/measure_nsc_atlas_dirac_preparation.py')


def digest(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def evaluate(cells,mass,region,k):
    with threadpool_limits(limits=1):
        start=time.process_time()
        chart=CompactMassChart(WallRadarChart(),mass)
        S=chart.duration
        # Parallelism is across eight independent cases, not nested pools.
        basis=ClockEnergyBasis(chart,cells,4,preparation_duration=S,
            region=region,angular_k=k,workers=1)
        if basis.flux_gram_defect>5e-12:
            raise RuntimeError('Characteristic transfer failed the joint flux Gram control')
        U,V=.8*S,.7*S
        radius=float(basis.points(U,V)['radius'])
        rows=[]
        matrices={}
        for count in (2,4):
            state=basis.select(count)
            if state['covariance'] is None:
                raise RuntimeError('Finite energy spectrum unresolved: '+str((cells,mass,region,k,count)))
            C=state['covariance']
            if np.linalg.norm(C@C-C,ord=2)>1e-10:
                raise RuntimeError('Finite covariance failed its purity/CAR control')
            rows.append(dict(radial_modes_per_leg=count,
                pressure_anisotropy=basis.pressure_difference(U,V,state),
                eigenvalues=state['eigenvalues'].tolist(),spectral_gap=state['gap'],
                relative_spectral_gap=state['gap']/float(np.linalg.norm(state['generator'],ord=2)),
                minimum_energy=state['minimum_energy'],
                covariance_purity_defect=float(np.linalg.norm(C@C-C,ord=2)),
                covariance_cross_leg_norm=float(np.linalg.norm(C[:count,count:],ord=2)),
                gram_condition=state['gram_condition'],gram_discrepancy=state['gram_discrepancy'],
                matrix_sha256=dict(generator=digest(state['generator']),covariance=digest(C))))
            matrices[count]=dict(generator=state['generator'],covariance=C)
        cpu=time.process_time()-start
        if cpu>30:
            raise RuntimeError('Case exhausted 30 CPU seconds')
        return dict(cells=cells,central_mass=mass,geometry=region,angular_k=k,
            duration=S,receiver_radius=radius,parent_mass=chart.parent_mass,
            flux_gram_defect=basis.flux_gram_defect,states=rows,_matrices=matrices,
            pid=os.getpid(),process_CPU_total=time.process_time())


def compare(rows):
    lookup={(r['cells'],r['central_mass'],r['geometry'],r['angular_k']):r for r in rows}
    controls=[];result={}
    for n in (64,128):
        for count in (2,4):
            def row(b,g,k):
                case=lookup[n,b,g,k]
                return next(s for s in case['states'] if s['radial_modes_per_leg']==count)
            differences={}
            for b in (0.,1e-4):
                differences[str(b)]=sum(row(b,'nested',k)['pressure_anisotropy']
                    -row(b,'cap',k)['pressure_anisotropy'] for k in (-1,1))
                for g in ('nested','cap'):
                    parity=abs(row(b,g,1)['pressure_anisotropy']-row(b,g,-1)['pressure_anisotropy'])
                    controls.append(dict(cells=n,radial_modes_per_leg=count,mass=b,geometry=g,
                        signed_sector_pressure_gap=parity))
            source_difference=differences['0.0001']-differences['0.0']
            radius=lookup[n,1e-4,'nested',1]['receiver_radius']
            A0=-3e-4/radius**3
            coefficient=8*np.pi*source_difference/A0
            result[f'{n}:{count}']=dict(
                nested_minus_cap_pressure_by_mass=differences,
                source_induced_pressure_difference=source_difference,
                source_induced_tidal_ratio_coefficient=coefficient,
                nominal_ratio_at_ell_10km=coefficient*constants.hbar*constants.G/constants.c**3/1e8,
                classical_tidal_denominator=A0)
    matrix_diagnostics=[]
    for b in (0.,1e-4):
        for g in ('nested','cap'):
            for k in (-1,1):
                for count in (2,4):
                    coarse=lookup[64,b,g,k]['_matrices'][count]
                    fine=lookup[128,b,g,k]['_matrices'][count]
                    state=next(s for s in lookup[128,b,g,k]['states']
                        if s['radial_modes_per_leg']==count)
                    difference=float(np.linalg.norm(fine['generator']-coarse['generator'],ord=2))
                    matrix_diagnostics.append(dict(mass=b,geometry=g,angular_k=k,
                        radial_modes_per_leg=count,operator_grid_spread=difference,
                        operator_spread_over_gap=difference/state['spectral_gap'],
                        covariance_grid_spread=float(np.linalg.norm(
                            fine['covariance']-coarse['covariance'],ord=2))))
    spreads=dict(tidal_grid_spread_by_mode_count={str(p):abs(
        result[f'128:{p}']['source_induced_tidal_ratio_coefficient']-
        result[f'64:{p}']['source_induced_tidal_ratio_coefficient']) for p in (2,4)},
        tidal_mode_count_spread_at_128=abs(result['128:4']['source_induced_tidal_ratio_coefficient']-
            result['128:2']['source_induced_tidal_ratio_coefficient']))
    return dict(finite_values=result,spreads=spreads,signed_sector_controls=controls,
        matrix_resolution=matrix_diagnostics,
        field_state_resolution=('resolved_against_two_grid_operator_spread' if all(
            r['operator_spread_over_gap']<.2 for r in matrix_diagnostics) else
            'unresolved_against_two_grid_operator_spread'),
        finite_algebra_controls='pass')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('unused compact output path with existing parent required')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    hashes={}
    for name in SOURCES:
        raw=(ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT)!=raw:
            raise ValueError('freeze the preparation method before execution: '+name)
        hashes[name]=hashlib.sha256(raw).hexdigest()
    started=time.perf_counter()
    rows=[];totals={};workers=min(8,os.cpu_count() or 1)
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        pending=[pool.submit(evaluate,n,b,g,k) for n in (64,128) for b in (0.,1e-4)
                 for g in ('nested','cap') for k in (-1,1)]
        try:
            for future in as_completed(pending):
                value=future.result()
                totals[value['pid']]=max(totals.get(value['pid'],0.),value['process_CPU_total'])
                if time.process_time()+sum(totals.values())+8>CPU_LIMIT:
                    raise RuntimeError('Shared preparation-reference CPU allowance exhausted')
                rows.append(value)
                print(json.dumps(dict(cells=value['cells'],mass=value['central_mass'],
                    geometry=value['geometry'],k=value['angular_k'])),flush=True)
        except Exception:
            for future in pending:
                future.cancel()
            raise
    comparison=compare(rows)
    for row in rows:
        row.pop('_matrices')
        row.pop('pid');row.pop('process_CPU_total')
    result=dict(schema='nsc-atlas-dirac-clock-energy-candidate-v1',method_commit=commit,
        sources_sha256=hashes,status='finite_candidate_only_UV_and_full_Psi_Phi_open',
        state_rule='negative projector of spacetime-windowed full clock energy; common peak-one bump',
        domain=dict(wall_clock_duration=rows[0]['duration'],preparation_duration=rows[0]['duration'],
            detector_U_over_S=.8,detector_V_over_S=.7,angular_sectors=[-1,1],
            magnetic_degeneracy_included_once=True,radial_modes_per_leg=[2,4],cells=[64,128],
            masses=[0.,1e-4],physical_reference_length_m=10000.,
            smooth_control='cap continuation in the same wall clock and compact-mass family'),
        cases=sorted(rows,key=lambda r:(r['cells'],r['central_mass'],r['geometry'],r['angular_k'])),
        comparison=comparison,
        limits=['No combined angular/frequency cutoff limit or local Hadamard proof.',
            'Point stress values are finite-candidate diagnostics, not renormalized vacuum stress.',
            'B/A is the local tidal channel, not the full Bardeen potential ratio Psi/Phi.',
            'The physical preparation rule is supplied, not derived from operator inheritance.',
            'P001 remains unconverged; no P001 evolution or new grid was run.'],
        compute=dict(workers=workers,BLAS_threads_per_worker=1,
            process_CPU_observed=time.process_time()+sum(totals.values()),
            process_CPU_charged=time.process_time()+sum(totals.values())+8,
            CPU_limit=CPU_LIMIT,wall_seconds=time.perf_counter()-started,
            array_memory_ceiling_bytes=2*1024**3,retained_full_arrays=False))
    raw=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>48000:
        raise RuntimeError('Compact preparation summary exceeds 48000 bytes; no data file saved')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),bytes=len(raw.encode()),
        compute=result['compute'],spreads=comparison['spreads'])),flush=True)


if __name__=='__main__':
    main()
