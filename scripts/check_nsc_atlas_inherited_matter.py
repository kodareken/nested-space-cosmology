"""Freeze rational inherited-background bounds, then evaluate the unaveraged target."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

from recursive_horizons.nsc_atlas_inherited_matter import (
    exact_initial_geometry, uniform_bound, dust_lightcone,
    scalar_lightcone_bound, outward_decimal_interval, solve_parent,
)

ROOT=Path(__file__).resolve().parents[1]
FORECAST=ROOT/'results/summaries/nsc-atlas-inherited-matter-forecast-v1.json'
VALIDATION=ROOT/'results/summaries/nsc-atlas-inherited-matter-validation-v1.json'
SOURCES=('src/recursive_horizons/nsc_atlas_inherited_matter.py',
    'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
    'src/recursive_horizons/nsc_atlas_flux_conversion.py',
    'scripts/check_nsc_atlas_inherited_matter.py')


def sources():
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def forecast():
    rows=[]
    for w in (20,2000,10**9):
        b=uniform_bound(w);light=scalar_lightcone_bound(w)
        lower_A=(F(4,3)*b['lower'])**F(1,3)
        upper_A=(F(4,3)*b['upper'])**F(1,3)
        rows.append(dict(omega=w,
            comoving_energy_factors_exact={k:str(b[k]) for k in ('delta','Gamma','lower','upper')},
            density_relative_error_ceiling=float(b['density_relative_error']),
            H_squared_relative_error_ceiling=float(b['hubble_squared_relative_error']),
            Lambda_DL_squared_relative_error_ceiling=float(light['relative_squared_distance_error']),
            Lambda_DL_squared_interval=outward_decimal_interval(light['lower'],light['upper']),
            energy_equality_scale_interval_approximate=[float(lower_A),float(upper_A)]))
    geometry=exact_initial_geometry();reference=dust_lightcone()
    record=dict(schema='nsc-atlas-inherited-matter-forecast-v1',
        result='One inherited scalar has a uniform pressureless-background comparison and a rational light-cone enclosure, without an independently supplied dust normalization.',
        exact_initial_geometry={k:str(v) for k,v in geometry.items()},
        ancestry=dict(field='uniform grandchild-minimum field; zero scalar momentum',
            fluxes=['1/2','1','3/2'],mixing='1/4',charge='1/2',tension='1/5',
            inherited_excitation='(1-epsilon)*(q_parent-q_grandchild)^2/2=3/8',
            inherited_vacuum='epsilon*q_parent^2/2=9/32',initial_excitation_over_vacuum='4/3'),
        theorem=dict(domain='A>=1 in the protected homogeneous closed parent, turning-point initial preparation',
            corrected_energy='I=A^3 E+(3/2)H A^3 chi Pi; I_dot=(3/2)H_dot A^3 chi Pi',
            constants=dict(H_lower='133/500',H_upper='117/250',H_initial_upper='443/1000',curvature_upper='57/2500'),
            rational_comparison='delta=3H_upper/(2omega); Gamma=3V_H/(2omega(1-delta)); L=(1-Gamma)/(1+delta); U=1/((1-Gamma)(1-delta))',
            pressure='Instantaneous pressure is retained and can oscillate with order-E amplitude; the theorem controls density, H and distances, not H_dot or the full perturbation/halo reduction.'),
        observer=dict(observer_scale=2,redshift='1/4',emitter_scale='8/5',
            motion='comoving negligible test source and receiver at parent pole',
            observable='Lambda_v D_L^2; dimensionless, c=1; no SI ruler',
            geometry='closed FLRW D_L=a0 A_observer (1+z) sin(eta)'),
        dust_comparison=dict(Lambda_DL_squared_interval=outward_decimal_interval(reference['lower'],reference['upper']),
            binomial_terms=reference['terms'],binomial_max_abs_argument=float(reference['binomial_q']),
            integral_tail_ceiling=float(reference['integral_tail']),
            method='Exact rational polynomial integration, geometric-series tail, alternating sinc bounds; strings are authoritative outward-rounded decimal intervals.'),
        frequency_members=rows,
        calibration='none; existing declared action/state inputs precede all target evaluations',
        measured_observations_used=False,absolute_vacuum_selected=False,
        independent_unaveraged_target_evaluated=False,
        source_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_sha256=sources())
    if FORECAST.exists():
        raise RuntimeError('creation-only forecast already exists')
    FORECAST.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(forecast=str(FORECAST.relative_to(ROOT)),dust=record['dust_comparison'],members=rows),indent=2))


def validate():
    if VALIDATION.exists():raise RuntimeError('creation-only validation already exists')
    frozen=json.loads(FORECAST.read_text(encoding='utf-8'))
    # The analytic predictor must still produce the frozen exact interval.
    for row in frozen['frequency_members']:
        light=scalar_lightcone_bound(row['omega'])
        assert row['Lambda_DL_squared_interval']==outward_decimal_interval(light['lower'],light['upper'])
    jobs=[dict(clock='proper',omega=20,rtol=2e-11),
          dict(clock='efolds',omega=20,rtol=2e-12),
          dict(clock='proper',omega=2000,rtol=2e-11,short_cycles=4)]
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(solve_parent,**job) for job in jobs]
        rows=[f.result() for f in futures]
    calls=sum(r['rhs_calls'] for r in rows)
    cpu=sum(r['process_cpu_seconds'] for r in rows)
    if calls>180000 or cpu>180:raise RuntimeError('declared calculation allowance exhausted')
    if max(r['maximum_Friedmann_residual'] for r in rows)>3e-9:
        raise RuntimeError('unaveraged gravitational constraint failed')
    if max(r['maximum_energy_work_residual'] for r in rows)>3e-9:
        raise RuntimeError('unaveraged scalar work balance failed')
    gaps={k:abs(rows[0][k]-rows[1][k]) for k in
          ('Lambda_DL_squared','energy_equality_scale','final_comoving_energy_ratio')}
    if max(gaps.values())>2e-8:raise RuntimeError('independent clocks disagree')
    record=dict(schema='nsc-atlas-inherited-matter-validation-v1',
        forecast_path=FORECAST.relative_to(ROOT).as_posix(),
        forecast_sha256=hashlib.sha256(FORECAST.read_bytes()).hexdigest(),
        unaveraged_rows=rows,independent_clock_gaps=gaps,
        total_rhs_calls=calls,total_solver_process_cpu_seconds=cpu,
        long_frequency_1e9_evolution_run=False,
        physical_observation_test=False,
        source_sha256=sources())
    VALIDATION.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--validate',action='store_true')
    args=p.parse_args()
    validate() if args.validate else forecast()
