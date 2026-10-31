#!/usr/bin/env python3
"""Compact regeneration of the same-action one-mode initial-neck result."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_circulation import Parent
from recursive_horizons.nsc_atlas_valve_neck import solve_necks,neck_diagnostics,critical_amplitude


def calculation():
    started=time.perf_counter();p=Parent.retained();z=p.mass_squared/p.H**2
    with threadpool_limits(limits=1):
        levels=[neck_diagnostics(data) for L in (8,16,24,32) for data in solve_necks(L,parent=p)]
        critical=[critical_amplitude(L,parent=p) for L in (16,24,32)]
    root=Path(__file__).resolve().parents[1]
    files=('src/recursive_horizons/nsc_atlas_valve_neck.py','scripts/check_nsc_atlas_valve_neck.py',
           'src/recursive_horizons/nsc_atlas_parent_circulation.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py')
    finest=[row for row in levels if row['degree']==32 and row['epsilon']==1.][0]
    return dict(schema='nsc-atlas-valve-neck-v1',
        model=dict(action='Unchanged retained canonical volume scalar plus real massive scalar',
            preparation='chi=epsilon*A*X1, s=sstar, all scalar normal velocities=0, K_ij=0',
            changed_input='All round-reference excitation energy placed in the lowest odd mode; previous equal two-mode occupation was a supplied preparation.',
            H_reference=p.H,vacuum=p.vacuum,mass_squared=p.mass_squared,kappa=p.kappa,
            z=z,kappa_A_squared=24/(z+3)),
        identities=dict(constraint='-8 Delta0 psi+6psi-kappa eps^2 A^2(1-u^2)psi-(6+kappa eps^2 A^2 z u^2)psi^5=0',
            neck='u=X1=0, reflection-fixed totally geodesic S2 in h=psi^4 H^(-2)gamma',
            slope_direct='dot(theta_pm)/H^2=psi_neck^(-4)*(2-4psi_uu/psi-kappa eps^2 A^2)',
            slope_reduced='dot(theta_pm)/H^2=3-(1+kappa eps^2 A^2/2)*psi_neck^(-4)',
            local_development='For an exact nearby smooth constraint solution with the negative sign retained, symmetry and local development imply a strictly future-trapped neck at sufficiently small positive proper time.'),
        levels=levels,critical_amplitudes=critical,
        convergence=dict(critical_amplitude_adjacent_spreads=[abs(critical[i+1]['epsilon_crit']-critical[i]['epsilon_crit']) for i in range(2)],
            critical_full_conformal_RMS_adjacent_ratios=[critical[i]['full_conformal_PDE_RMS']/critical[i+1]['full_conformal_PDE_RMS'] for i in range(2)],
            epsilon1_neck_slope_L24_L32_spread=abs(next(row for row in levels if row['degree']==24 and row['epsilon']==1.)['dot_theta_pm_over_H_squared_direct']-finest['dot_theta_pm_over_H_squared_direct'])),
        scope=dict(numerically_resolved_initial_constraint=True,certified_continuum_enclosure=False,
            spacetime_evolution_performed=False,regular_child_proved=False,received_child_energy_calculated=False,
            occupation_selected_by_inheritance=False,SI_scale_selected=False),
        reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_neck.py --output results/summaries/nsc-atlas-valve-neck.json',
        source_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files},
        compute=dict(BLAS_threads=1,elapsed_seconds=time.perf_counter()-started,raw_arrays_saved=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('results/summaries/nsc-atlas-valve-neck.json'))
    args=parser.parse_args();result=calculation()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(args.output),critical_amplitudes=[dict(degree=r['degree'],epsilon_crit=r['epsilon_crit']) for r in result['critical_amplitudes']],elapsed_seconds=result['compute']['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
