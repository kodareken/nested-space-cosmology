"""Conserved local scalar structure and exact exterior matching."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from recursive_horizons.nsc_atlas_core_reservoir import CoreReservoir

ROOT=Path(__file__).resolve().parents[1]


def identities():
    R,k,U0,U,M,Mp,v,vp,Pi,F,b=sp.symbols('R k U0 U M Mp v vp Pi F b',real=True)
    B=1+v*v-k*U0*R*R/3-2*M/R
    Bp=2*v*vp-2*k*U0*R/3-2*Mp/R+2*M/R**2
    Kr=-vp-k*R*Pi*F/2;Kt=-v/R
    rho=(Pi*Pi+B*F*F)/2+U
    momentum=-2*b*(-vp/R+v/R**2)+2*b*(Kr-Kt)/R+k*Pi*b*F
    ham=-2*Bp/R+2*(1-B)/R**2+4*Kr*Kt+2*Kt*Kt-2*k*rho
    # Differentiate the geometric definition using the radial Einstein
    # acceleration and the mixed Einstein equation for n(beta).
    m=sp.symbols('m',real=True)
    pr=(Pi*Pi+b*b*F*F)/2-U
    nv=-m/R**2-k*R*pr/2;nb=-k*R*b*Pi*F/2
    nm=v*m/R+R*(v*nv-b*nb)
    target=k*R*R*(Pi*b*b*F-v*pr)/2
    mr=k*R*R*((Pi*Pi+b*b*F*F)/2+U-v*Pi*F)/2
    fixed=target-v*mr
    raw=dict(momentum_constraint=momentum,
             Hamilton_mass_constraint=ham.subs(Mp,k*R*R*(rho-U0-v*Pi*F)/2),
             normal_mass_transport=nm-target,
             outgoing_fixed_radius_flux=sp.expand(fixed.subs(Pi,-b*F)+k*R*R*b*F*F*(b+v)**2/2),
             homogeneous_pressure_work=target.subs({F:0,Pi:0})-k*R*R*v*U/2,
             fixed_radius_vacuum_work_cancels=fixed.subs({F:0,Pi:0}))
    ans={name:str(sp.simplify(value)) for name,value in raw.items()}
    if any(v!='0' for v in ans.values()):raise RuntimeError(ans)
    return ans


def worker(label):
    start=time.process_time()
    config=CoreReservoir(core_enabled=label!='zero',pulse_fraction=0. if label in ('zero','no_pulse') else .02)
    sol=config.solve(rtol=2e-11 if label=='refined' else 2e-10,
                     substeps=240 if label=='refined' else 150)
    summary=sol.summary()
    radius=np.linspace(0,config.outer,2001)
    state=sol.state(radius)
    columns=np.array([radius,sol.solution.sol(radius)[0],state['field_displacement'],
                      state['Pi'],state['beta']]).T
    summary['sampled_state_sha256']=hashlib.sha256(np.asarray(columns,dtype='<f8').tobytes()).hexdigest()
    summary['cpu_seconds']=time.process_time()-start
    return label,summary,columns


def calculation():
    start=time.process_time();wall=time.perf_counter()
    symbolic=identities()
    with ProcessPoolExecutor(max_workers=4) as pool:
        results=list(pool.map(worker,('compensated','refined','no_pulse','zero')))
    runs={name:r for name,r,_ in results};main=runs['compensated'];fine=runs['refined']
    baseline=json.loads((ROOT/'results/summaries/nsc-atlas-scalar-core-v1.json').read_text(encoding='utf-8'))
    if not np.isclose(main['weak_amplitude'],baseline['physical_hierarchy']['formed_coherent_weak_amplitude'],atol=1e-14,rtol=0):
        raise RuntimeError('core representative no longer matches the cooling record')
    comparison={k:fine[k]-main[k] for k in ('depletion_fraction','core_geometric_excess',
        'pulse_geometric_excess','proper_energy_excess','parent_proper_length_change')}
    files=('src/recursive_horizons/nsc_atlas_core_reservoir.py','src/recursive_horizons/nsc_atlas_scalar_core.py',
           'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py','src/recursive_horizons/nsc_atlas_flux_conversion.py',
           'src/recursive_horizons/nsc_atlas_flux_transition.py','scripts/check_nsc_atlas_core_reservoir.py',
           'scripts/plot_nsc_atlas_core_reservoir.py','tests/test_nsc_atlas_core_reservoir.py',
           'results/summaries/nsc-atlas-scalar-core-v1.json')
    cpu=time.process_time()-start+sum(r['cpu_seconds'] for r in runs.values())
    record=dict(schema='nsc-atlas-core-reservoir-v1',
        result='A nonlinear local core, inherited-scalar depletion and outward-current pulse satisfy the shared atlas constraints with zero enclosing mass offset and unchanged initial wall data.',
        symbolic_zero_identities=symbolic,runs=runs,independent_tolerance_comparison=comparison,
        preparation=dict(core='Same mass-selected coherent representative as scalar-core-v1; not its full cooling snapshot.',
            depletion='-a*chi0*1024*x^5*(1-x)^5 on [Rc,2Rc], a solved by M(4Rc)=0.',
            pulse='.02*chi0*1024*x^5*(1-x)^5*sin(omega*(R-3Rc)) on [3Rc,4Rc]; Pi=-beta*phi_R.',
            slicing='v=H_reference*R; Kt=-H_reference; Kr=-H_reference-kappa*R*Pi*phi_R/2.',
            supplied='Frequency, core norm/phase, radial supports and pulse shape/amplitude are preparation inputs. Only depletion is fixed by compensation.'),
        equations=dict(mass='m=R*(1+v^2-beta^2)/2=M+kappa*U0*R^3/6',
            radial='M_R=kappa*R^2*(rho-U0-v*Pi*phi_R)/2',
            normal='n(m)=kappa*R^2*(Pi*beta^2*phi_R-v*p_r)/2',
            fixed_radius_outgoing='(n-v*s/beta)(m)=-kappa*R^2*beta*phi_R^2*(beta+v)^2/2 for Pi=-beta*phi_R'),
        limits='Initial nonlinear compatibility plus exact transport identities; no finite formation evolution, unique seed, abundance, nonspherical stability or observational fit. Exterior initial data coincide; equality of future observables is restricted to their common causal domain. The old homogeneous-parent shielding bound is not assigned to the altered interior.',
        next='Evolve the compensated scalar state and geometry together, retaining the outgoing field and pressure work; derive the initial cloud from inheritance and connect to an observer signal.',
        worker_count=4,process_cpu_seconds_including_workers=cpu,wall_seconds=time.perf_counter()-wall,
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_core_reservoir.py --figure')
    if cpu>120:raise RuntimeError('120 CPU-second compensation budget exceeded')
    for row in runs.values():
        if max(abs(row['total_geometric_excess']),abs(row['independent_mass_gap']))>2e-11:
            raise RuntimeError('independent compensation/mass check failed')
        if row['maximum_sampled_Hamilton_residual']>1e-7 or row['maximum_sampled_momentum_residual']>1e-12:
            raise RuntimeError('Einstein constraints failed')
    if abs(comparison['depletion_fraction'])>1e-7:raise RuntimeError('depletion tolerance unresolved')
    return record,results[0][2]


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path);parser.add_argument('--figure',action='store_true');args=parser.parse_args()
    result,columns=calculation();encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>20000:raise RuntimeError('20 KB compact record limit exceeded')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
    if args.figure:
        from plot_nsc_atlas_core_reservoir import render
        render(result,columns)
