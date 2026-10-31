"""Reproduce finite future-trapped-neck development and the regional ledger."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve

ROOT=Path(__file__).resolve().parents[1]


def band_ledger(system,solution,t=.2,step=2e-5):
    def state(at):return solution.sol(at).reshape(8,system.modes)
    band=system.trapped_band(state(t))
    if band is None:return None
    earlier=system.trapped_band(state(t-step));later=system.trapped_band(state(t+step))
    du=(later['u_edge']-earlier['u_edge'])/(2*step)
    derivative=(later['energy']-earlier['energy'])/(2*step)
    area=4*math.pi*band['boundary_radius']**2
    speed=band['boundary_radial_scale']*du/math.sqrt(1-band['u_edge']**2)
    physical_influx=-2*area*band['outward_energy_flux']
    sweep=2*area*band['boundary_rho']*speed
    predicted=band['pressure_work_rate']+physical_influx+sweep
    return dict(proper_time=t,finite_difference_step=step,trapped_band=band,
                energy_derivative=derivative,pressure_work_rate=band['pressure_work_rate'],
                material_influx_rate=physical_influx,boundary_sweep_rate=sweep,
                boundary_speed_relative_to_normal=speed,
                sum_of_ledger_terms=predicted,ledger_absolute_residual=abs(derivative-predicted),
                meaning='A moving marginal boundary can enclose pre-existing energy. Sweep is not extra parent-to-child material transfer.')


def calculation():
    started=time.perf_counter();levels=[]
    sources=('src/recursive_horizons/nsc_atlas_valve_neck.py','src/recursive_horizons/nsc_atlas_valve_spherical.py',
             'src/recursive_horizons/nsc_atlas_parent_circulation.py','src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
             'scripts/check_nsc_atlas_valve_spherical.py')
    sha=lambda name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    before={name:sha(name) for name in sources}
    with threadpool_limits(limits=1):
        for modes in (24,32,48):
            system=SphericalValve(modes);solution=system.evolve(.3)
            def state(t):return solution.sol(t).reshape(8,modes)
            snapshots=[system.diagnostic(t,state(t)) for t in (0.,.02,.1,.3)]
            bands=[dict(proper_time=t,**system.trapped_band(state(t))) for t in (.02,.1,.3)]
            x,w=leggauss(32)
            work=.15*sum(float(weight)*system.diagnostic(float(t),state(t))['pressure_work_rate'] for t,weight in zip(.15*(x+1),w))
            delta=snapshots[-1]['energy']-snapshots[0]['energy']
            levels.append(dict(modes=modes,snapshots=snapshots,trapped_bands=bands,
                full_region_work=dict(energy_change=delta,independent_integrated_pressure_work=work,
                                      residual=abs(delta-work),half_region_received_material_energy=0.,
                                      reason='Reflection makes the total physical timelike energy flux at X1=0 identically zero.'),
                moving_band_ledger=band_ledger(system,solution),rhs_calls=solution.nfev))
        fine=levels[-1]
        comparison=[dict(proper_time=t,
            theta_pm_over_H=[level['snapshots'][i]['theta_plus_over_H'] for level in levels],
            adjacent_spreads=[abs(levels[j+1]['snapshots'][i]['theta_plus_over_H']-levels[j]['snapshots'][i]['theta_plus_over_H']) for j in (0,1)],
            Hamiltonian_RMS_over_H_squared=[level['snapshots'][i]['Hamiltonian_RMS_over_H_squared'] for level in levels])
            for i,t in enumerate((0.,.02,.1,.3))]
    return dict(schema='nsc-atlas-valve-spherical-v1',source_sha256=before,
        producer_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_unchanged=before=={name:sha(name) for name in sources},
        preparation='Declared new one-mode chi=epsilon*A*X1, epsilon=1, same round-reference energy normalization and action; not selected by inheritance.',
        equations=dict(metric='ds²=-dt²+q² dchi²+b² sin²chi dOmega²; u=coschi; q=b+(1-u²)c',
            curvature='K_rad=K_tan+(1-u²)d',
            future_null_expansions='theta_pm=-2 K_tan +/- 2 R_chi/(q R)',
            regional_ledger='dE/dt=pressure work + physical material influx + moving-boundary sweep',
            mass='m_geo=R/2*(1+K_tan² R²-(R_chi/q)²)',
            interior_frame='theta_plus*theta_minus>0 implies grad(R) timelike. -R is a local future time inside the strict trap; synchronous variables remain regular through the marginal boundary.'),
        numerical=dict(method='DOP853',proper_endpoint=.3,rtol=3e-11,atol=1e-13,max_step=.002,
            initial_degree=40,spatial_modes=[24,32,48],BLAS_threads=1,elapsed_seconds=time.perf_counter()-started,
            raw_arrays_stored=False),levels=levels,convergence=comparison,
        verdict=dict(resolved_future_trapped_neck=True,global_event_horizon_proved=False,
            earlier_formation_history_selected=False,rotation_used=False,external_compression_supplied=False,
            internal_preparation_changes=True,regular_child_proved=False,one_way_child_material_supply_proved=False,
            child_vacuum_and_equation_of_state=None,child_energy_comparable_to_Parke_three_eighths=None,
            old_8847_percent_applies=False,
            next='Follow the actual receiving trapped band and its horizons into a regular continuation. Zero center flux does not exclude symmetric horizon influx. Retain inherited state; do not repeat this finite threshold calculation.'),
        reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_spherical.py --output results/summaries/nsc-atlas-valve-spherical-v1.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('results/summaries/nsc-atlas-valve-spherical-v1.json'))
    args=parser.parse_args();result=calculation();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(args.output),convergence=result['convergence'],ledger=result['levels'][-1]['moving_band_ledger']),indent=2))


if __name__=='__main__':main()
